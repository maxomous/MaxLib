#include "String.h"
    
using namespace std;
    
namespace MaxLib {
namespace String {
                
    string va_str(const char* format, ... )
    { 
        va_list arglist;
        char buf[255];
        va_start( arglist, format );
        vsnprintf(buf, sizeof(buf), format, arglist);
        va_end( arglist );
        
        return string(buf);
    }
      
    void LowerCase(string& str) {
        transform(str.begin(), str.end(),str.begin(), ::tolower);
    }

    void UpperCase(string& str) {
        transform(str.begin(), str.end(),str.begin(), ::toupper);
    }
    
    
} // end namespace String
} // end namespace Maxomous
