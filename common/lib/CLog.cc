#include "CLog.hh"
#include <cstdarg>
#include <cstdio>

namespace SiFi {
namespace log {

void info(const char* fmt, ...) {
  va_list args;
  printf(fmt, args);
}
void debug(const char* fmt, ...) {
  va_list args;
  printf(fmt, args);
}
void warn(const char* fmt, ...) {
  va_list args;
  printf(fmt, args);
}
void error(const char* fmt, ...) {
  va_list args;
  printf(fmt, args);
}

} // namespace log
} // namespace SiFi
