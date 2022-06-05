/*
 * file.hpp
 */

#pragma once

#include <string>
#include <vector>
#include <functional>
#include <mutex>


namespace MaxLib {
namespace File {

struct FileDesc {
    enum Type { None, Folder, File };
    Type type;                      // file or folder
    std::string name;               // filename / foldername
    std::string ext;                // file extension / "" for folder
    tm lastModified;                // last modification
    std::string lastModifiedAsText; // last modified as a readable string
    int id;                         // ID to allow us to sort a list of FileDesc's and find the right item
};


// Returns the directory of this executable
std::string ThisDir();
// Returns the directory of this executable and appends location to it
std::string ThisDir(const std::string& location);
// Combines dir & name to give a file location
std::string CombineDirPath(const std::string& dir, const std::string& name);
// Opens the file dialog
int OpenFileDialog();
// checks if file already exists
bool Exists(const std::string& name);
// Writes to a new file or overwrites existing file
void Write(const std::string& filename, const std::string& str);
// Writes to a new file or appends to the end of an existing file
void Append(const std::string& filename, const std::string& str);
// Writes many strings in one go to a file
void WriteArray(const std::string& filename, const std::vector<std::string>& array);
// Reads contents of file and calls func callback
// returns 0 on success / 1 if unsuccessful
// Usage:
//    auto executeLine = [](string& str) {
//        cout << str << endl;
//        return 0;
//    }; 
//    if(File::Read("/home/pi/Desktop/New.nc", executeLine)) {
//        cout << "Error: Could not open file" << endl;
//    }

int Read(const std::string& filename, const std::function<int(std::string&)>& callback, uint firstLine = 1, uint lastLine = 0);
// returns the number of lines in a file
// returns -1 on failure
int GetNumLines(const std::string& filename);
// retrieves information about the files within a given directory
// if extensions is not an empty string, it will only return files with given extensions (can be seperated by ',' e.g. "exe,ini")
// returns 1 on failure
int GetFilesInDir(const std::string& location, const std::string& extensions, std::vector<FileDesc>& files);

} // end namesapce File
} // end namesapce MaxLib
