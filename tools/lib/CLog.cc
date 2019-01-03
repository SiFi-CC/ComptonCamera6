#include "CLog.hh"
#include <cstdarg>
#include <cstdio>
#include <string>

namespace SiFi {
namespace log {

level currentLogLevel = level::info;

void setLevel(level logLevel) { currentLogLevel = logLevel; }
level getLevel() { return currentLogLevel; };

void debug(const char* fmt, ...) {
  if (currentLogLevel > level::debug) { return; }
  va_list args;
  va_start(args, fmt);
  vprintf((std::string(fmt) + "\n").data(), args);
  va_end(args);
}
void info(const char* fmt, ...) {
  if (currentLogLevel > level::info) { return; }
  va_list args;
  va_start(args, fmt);
  vprintf((std::string(fmt) + "\n").data(), args);
  va_end(args);
}
void warn(const char* fmt, ...) {
  if (currentLogLevel > level::warn) { return; }
  va_list args;
  va_start(args, fmt);
  vprintf((std::string(fmt) + "\n").data(), args);
  va_end(args);
}

void error(const char* fmt, ...) {
  va_list args;
  va_start(args, fmt);
  vprintf((std::string(fmt) + "\n").data(), args);
  va_end(args);
}

} // namespace log
} // namespace SiFi
