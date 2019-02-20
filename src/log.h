/*
 * log.h
 *
 *  Created on: May 26, 2014
 *      Author: meiermark
 */

#ifndef LOG_H_
#define LOG_H_

#include <sstream>
#include <string>
#include <cstdio>
#include <sys/time.h>

inline std::string NowTime();

enum LogLevel {ERROR, WARNING, INFO, DEBUG, DEBUG1, DEBUG2, DEBUG3, DEBUG4};

class Log
{
public:
    Log();
    virtual ~Log();
    std::ostringstream& Get(LogLevel level = INFO);
public:
    static LogLevel& reporting_level();
    static std::string to_string(LogLevel level);
    static LogLevel from_string(const std::string& level);
    static LogLevel from_int(const int level);
protected:
    std::ostringstream os;
private:
    Log(const Log&);
    Log& operator =(const Log&);
};

inline Log::Log()
{}

inline std::ostringstream& Log::Get(LogLevel level)
{
    os << "- " << NowTime();
    os << " " << to_string(level) << ": ";
    os << std::string(level > DEBUG ? level - DEBUG : 0, '\t');
    return os;
}

inline Log::~Log()
{
    os << std::endl;
    fprintf(stderr, "%s", os.str().c_str());
    fflush(stderr);
}

inline LogLevel& Log::reporting_level()
{
    static LogLevel reportingLevel = DEBUG4;
    return reportingLevel;
}

inline std::string Log::to_string(LogLevel level)
{
	static const char* const buffer[] = {"ERROR", "WARNING", "INFO", "DEBUG", "DEBUG1", "DEBUG2", "DEBUG3", "DEBUG4"};
    return buffer[level];
}

inline LogLevel Log::from_string(const std::string& level)
{
    if (level == "DEBUG4")
        return DEBUG4;
    if (level == "DEBUG3")
        return DEBUG3;
    if (level == "DEBUG2")
        return DEBUG2;
    if (level == "DEBUG1")
        return DEBUG1;
    if (level == "DEBUG")
        return DEBUG;
    if (level == "INFO")
        return INFO;
    if (level == "WARNING")
        return WARNING;
    if (level == "ERROR")
        return ERROR;
    Log().Get(WARNING) << "Unknown logging level '" << level << "'. Using INFO level as default.";
    return INFO;
}

inline LogLevel Log::from_int(const int level) {
	  if (level == 7)
	    return DEBUG4;
	  if (level == 6)
	    return DEBUG3;
	  if (level == 5)
	    return DEBUG2;
	  if (level == 4)
	    return DEBUG1;
	  if (level == 3)
	    return DEBUG;
	  if (level == 2)
	    return INFO;
	  if (level == 1)
	    return WARNING;
	  if (level == 0)
	    return ERROR;
	  Log().Get(WARNING) << "Unknown logging level '" << level
	                     << "'. Using INFO level as default.";
	  return INFO;
}

typedef Log FILELog;

#define HH_LOG(level) \
  if (level <= Log::reporting_level()) \
    Log().Get(level)

inline std::string NowTime()
{
    char buffer[11];
    time_t t;
    time(&t);
    tm r = {0};
    strftime(buffer, sizeof(buffer), "%X", localtime_r(&t, &r));
    struct timeval tv;
    gettimeofday(&tv, 0);
    char result[100] = {0};
    sprintf(result, "%s.%03ld", buffer, (long)tv.tv_usec / 1000);
    return result;
}

#endif /* LOG_H_ */
