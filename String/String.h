#pragma once
#include <stdarg.h>     // va_list, va_start, va_arg, va_end
#include <algorithm>
#include <string>


namespace MaxLib {
namespace String {

   // converts variable arguments to a string
    std::string va_str(const char *format, ...);
    // modifies string to lower case
    void LowerCase(std::string &str);
    // modifies string to upper case
    void UpperCase(std::string &str);
    
} // end namespace String
} // end namespace MaxLib
