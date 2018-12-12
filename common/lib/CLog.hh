#ifndef __CLog_H_
#define __CLog_H_ 1

namespace SiFi {
namespace log {

void info(const char* fmt, ...);
void debug(const char* fmt, ...);
void warn(const char* fmt, ...);
void error(const char* fmt, ...);

} // namespace log
} // namespace SiFi

#endif
