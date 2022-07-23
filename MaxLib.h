#pragma once
/*
 * Name: 
 *    MaxLib
 * Description: 
 *    A Set of Utility functions (written on the Raspberry Pi 3b / 4)
 */

#include "MaxLib/ads1115.h"
#include "MaxLib/File.h"
#include "MaxLib/Geom.h"
#include "MaxLib/Vector.h"
#include "MaxLib/String.h"

namespace MaxLib {
        
    
  
// Normalises seconds into hours, minutes & seconds
class Time {
public:
    Time(uint hours, uint minutes, uint seconds) : m_hr(hours), m_min(minutes), m_sec(seconds) {}
    Time(uint seconds) { 
        m_hr = seconds / 3600;
        seconds %= 3600;
        m_min = seconds / 60;
        seconds %= 60;
        m_sec = seconds; 
    }
    uint Hours() { return m_hr; }
    uint Mins() { return m_min; }
    uint Secs() { return m_sec; }
    std::string TimeString() { return MaxLib::String::va_str("%u:%.2u:%.2u", m_hr, m_min, m_sec); }
private:
    uint m_hr;
    uint m_min;
    uint m_sec;
};



class Log {
public:
    enum LogLevel {
        LevelInfo,
        LevelWarning,
        LevelError,
        LevelCritical,
        LevelDebug
    };
    
    static std::string LevelPrefix(LogLevel level) {
        if (level == LevelInfo)
            return "[Info] ";
        else if (level == LevelWarning)
            return "[Warning] ";
        else if (level == LevelError)
            return "[Error] ";
        else if (level == LevelCritical)
            return "[Critical] ";
        else if (level == LevelDebug)
            return "[Debug] ";
        return "";
    }
    
    static void SetLevel(LogLevel level) {
        // lock the mutex
        std::lock_guard<std::mutex> guard(get().m_mutex);
        get().m_logLevelTerminal = level;
    }
    // set bit flags for use with debugging
    static void SetDebugFlags(int flags) {
        // lock the mutex
        std::lock_guard<std::mutex> guard(get().m_mutex);
        get().m_debugFlags = flags;
    }
    // returns unique id of handler
    static void RegisterHandler(const std::function<void(const char*, LogLevel, const std::string&)>& eventHandler) {
        std::lock_guard<std::mutex> guard(m_mutex);
        get().m_PrintHandler = eventHandler;
    }
    static void Debug(int flag, const std::string &msg) { get().Print(flag, LevelDebug, msg.c_str()); }
    static void Critical(const std::string &msg) { get().Print(LevelCritical, msg.c_str()); }
    static void Error(const std::string &msg) { get().Print(LevelError, msg.c_str()); }
    static void Warning(const std::string &msg) { get().Print(LevelWarning, msg.c_str()); }
    static void Info(const std::string &msg) { get().Print(LevelInfo, msg.c_str()); }

    template <typename... Args>
    static void Debug(int flag, const char *msg, Args... args) { get().Print(flag, LevelDebug, msg, args...); }
    template <typename... Args>
    static void Critical(const char *msg, Args... args) { get().Print(LevelCritical, msg, args...); }
    template <typename... Args>
    static void Error(const char *msg, Args... args) { get().Print(LevelError, msg, args...); }
    template <typename... Args>
    static void Warning(const char *msg, Args... args) { get().Print(LevelWarning, msg, args...); }
    template <typename... Args>
    static void Info(const char *msg, Args... args) { get().Print(LevelInfo, msg, args...); }

private:

    LogLevel m_logLevelTerminal = LevelInfo; // default show all
    int m_debugFlags = 0;                    // default show none
    std::mutex m_mutex;
    // user settable print handler
    std::function<void(const char*, LogLevel, const std::string&)> m_PrintHandler;

    
    template <typename... Args>
    void Print(int debugFlag, LogLevel level, const char *msg, Args... args) {
        if (!PrintDebug(debugFlag))
            return;
        Print(level, msg, args...);
    }


    template <typename... Args>
    void Print(LogLevel level, const char *msg, Args... args) {
        // lock the mutex
        std::lock_guard<std::mutex> guard(m_mutex);
        
        // Make date string
        char date[32];
        time_t t = time(NULL);
        strftime(date, 32, "[%H:%M:%S]", localtime(&t));
        // Make message string
        std::string message = MaxLib::String::va_str(msg, args...);
        
        // Print to terminal
        if(m_logLevelTerminal <= level) { printf("%s%s%s\n", date, LevelPrefix(level).c_str(), message.c_str()); }
        // print to handler
        if(m_PrintHandler) { m_PrintHandler(date, level, message); }
        // stop program execution
        if (level == LevelCritical)
            exit(1);
    }


    bool PrintDebug(int flag) {
        std::lock_guard<std::mutex> guard(m_mutex);
        return m_debugFlags & flag;
    }

    
    static Log &get() {
        static Log log;
        return log;
    }

    Log() {}                                // delete the constructor
    Log(const Log &) = delete;              // delete the copy constructor
    Log &operator=(const Log &) = delete;   // delete the copy assignment operatory
};

} // end namespace MaxLib
