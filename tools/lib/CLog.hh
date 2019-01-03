#ifndef __CLog_H_
#define __CLog_H_ 1

namespace SiFi {
namespace log {

enum class level : int { debug = 1, info = 2, warn = 3, error = 4 };

void setLevel(level logLevel);
level getLevel();

/**
 * Loggs very verbose informations, data for every iteration, used to debug
 * problems with code.
 */
void debug(const char* fmt, ...);

/**
 * Logs basic information about current program, amount of this logs should not
 * clutter the screen
 */
void info(const char* fmt, ...);

/**
 * Log warning, should be used when sth unexpected happens but it's not critical
 * to program workings(e.g. particle with negative energy)
 */
void warn(const char* fmt, ...);

/**
 * Log errors, should be used if problem is critical and only solution is to
 * exit the program.
 */
void error(const char* fmt, ...);

} // namespace log
} // namespace SiFi

#endif
