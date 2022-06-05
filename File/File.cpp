
#include <cstring>
#include <algorithm>
#include <functional>
#include <iostream>
#include <fstream>
#include <sstream>
#include "dirent.h"  
#include <sys/stat.h> 
//#include <limits.h>
#include <unistd.h>

#include "File.h"

using namespace std;

namespace MaxLib {
namespace File {
    
// Mutext global to File
std::mutex m_mutex;

    
int OpenFileDialog() 
{
    string zenityCall = "/usr/bin/zenity";
    zenityCall += "  --file-selection --modal --title=\"Select File\"";
    
    FILE *f = popen(zenityCall.c_str(), "r");
    char buffer[1024];
    fgets(buffer, 1024, f);
    
    int ret = pclose(f);
    if(ret < 0)
        cerr << "Error: Couldn't open file" << endl;
    
    cout << "File Opened: " << buffer << endl;
    
    return ret;
}


string ThisDir() 
{
    // lock the mutex
    std::unique_lock<std::mutex> locker(m_mutex);
    char result[255];
    ssize_t count = readlink("/proc/self/exe", result, 255);
    string out = string(result, (count > 0) ? count : 0);
    return out.substr(0, out.find_last_of("/") + 1);
}

string ThisDir(const std::string& location)
{
    return ThisDir() + location;
}
    
string CombineDirPath(const string& dir, const string& name) {
    
    string filePath = dir;
    // append '/' if needed
    if(filePath.empty())        
        filePath += '/';
    if(filePath.back() != '/')  
        filePath += '/';
    filePath += name;
    return filePath;
}

bool Exists(const std::string& name) 
{
    ifstream f(name.c_str());
    return f.good();
}

void Write(const string& filename, const string& str) {
    // lock the mutex
    std::unique_lock<std::mutex> locker(m_mutex);
    
    ofstream f(filename, ios::out);
    if (f.is_open()) {
        f << str;
        f.close();
    } else {
        std::cout << "Error: " << filename << " is not open" << endl;
    }
}

void Append(const string& filename, const string& str) {
    // lock the mutex 
    std::unique_lock<std::mutex> locker(m_mutex);

    ofstream f(filename, ios::out | ios::app);
    if (f.is_open()) {
        f << str;
        f.close();
    } else {
        std::cout << "Error: " << filename << " is not open" << endl;
    }
}

void WriteArray(const string& filename, const vector<string>& array) {
    // lock the mutex 
    std::unique_lock<std::mutex> locker(m_mutex);

    ofstream f(filename, ios::out);
    if (f.is_open()) {
        for(const string& str : array) {
            f << str << '\n';
        }
        f.close();
    } else {
        std::cout << "Error: " << filename << " is not open" << endl;
    }
}

// reads a line of a file and invokes the callback function with a pointer to the string 
// returns 0 on success / 1 if unsuccessful
int Read(const string& filename, const function<int(string&)>& callback, uint firstLine, uint lastLine) 
{
    // lock the mutex
    std::unique_lock<std::mutex> locker(m_mutex);
    
    ifstream openFile(filename);
    if(!openFile.is_open()) {
        cout << "Error: Couldn't open file " << filename << endl;
        return -1;
    }
    string output;
    uint n = 0;
    while(getline(openFile, output)) {
        if(++n < firstLine)
            continue;
        if(lastLine != 0 && n > lastLine)
            break;
        if(callback(output)) {
            cout << "Error: Cannot execute line of file " << filename << endl;
            return -1;
        }
    }
    openFile.close();
    return 0;
}

int GetNumLines(const string& filename) 
{
    // lock the mutex
    std::unique_lock<std::mutex> locker(m_mutex);
    
    ifstream openFile(filename);
    if(!openFile.is_open()) {
        cout << "Error: Couldn't open file " << filename << endl;
        return -1;
    }
    string unused;
    int count = 0;
    while(getline(openFile, unused)) {
        count++;
    }
    openFile.close();
    return count;
}

// take a directory location and returns the files and directories within it in the vector files
// if extensions is not an empty string, it will only return files with given extensions (can be seperated by ',' e.g. "exe,ini")
// returns 1 on failure
int GetFilesInDir(const string& location, const string& extensions, vector<FileDesc>& files)
{
    // lock the mutex
    std::unique_lock<std::mutex> locker(m_mutex);
    
    auto buildExtensionsList = [](const string& ext) 
    {
        istringstream stream(ext);
        vector<string> list;
        // build vector of extensions
        for(string s; getline(stream, s, ','); ) {
            // strip out whitespace & add to vector
            s.erase(remove_if(s.begin(), s.end(), ::isspace), s.end());
            list.push_back(s);
        }
        return list;
    };
    
    auto isValidExtension = [](const string& ext, vector<string> list)
    {    // returns true if ext is part of list
        for(string s : list) {
            if (ext == s)
            return true;
        }
        return false;
    };
    
    files.clear();
    vector<string> permittedExt = buildExtensionsList(extensions);
    
    DIR *dir;
    struct dirent *ent;
    
    if ((dir = opendir (location.c_str())) != NULL) {
        while ((ent = readdir (dir)) != NULL) {
            // ignore . hidden files & ..
            if(strncmp(ent->d_name, ".", 1) && strncmp(ent->d_name, "..", 2)) 
            {
                static int id = 0;
                string str = ent->d_name;
                string name, ext;
                FileDesc::Type type = FileDesc::Type::None;
                // if directory
                if(ent->d_type == DT_DIR) {
                    name = str;
                    ext = "";
                    type = FileDesc::Type::Folder;
                } else {
                    // split name and file extension
                    size_t dotPos = str.find_last_of(".");
                    name = (dotPos != string::npos) ? str.substr(0, dotPos) : str;
                    ext = (dotPos != string::npos) ? str.substr(dotPos+1) : "";
                    type = FileDesc::Type::File;
                }
                if(type == FileDesc::Type::Folder || extensions == "" || isValidExtension(ext, permittedExt))
                {
                    struct stat attrib;
                    tm modtimeLocal;
                    char dateStr[32];
                    string filepath = CombineDirPath(location, string(ent->d_name));
                    // get modified date & time
                    if(!stat(filepath.c_str(), &attrib)) {
                        time_t modtime = attrib.st_mtime;
                        modtimeLocal = *localtime(&modtime);
                        strftime(dateStr, 32, "%d-%b-%y  %H:%M", localtime(&modtime));
                    }
                    files.push_back({ // FileDesc
                        type,         // .type
                        name,         // .name
                        ext,          // .ext
                        modtimeLocal, // .lastModified
                        dateStr,      // .lastModifiedAsText
                        id++,         // .id                            
                    });
                }
            }
        }

        closedir (dir);
    } else {
        /* could not open directory */
        perror ("");
        return 1;
    }
    return 0;
}

} // end namesapce File
} // end namesapce MaxLib
