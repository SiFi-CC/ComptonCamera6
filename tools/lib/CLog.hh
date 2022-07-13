#ifndef __CLog_H_
#define __CLog_H_ 1

#include <spdlog/spdlog.h>
#define __FORCE_SPDLOG_H_TO_BE_INCLUDED_FIRST_ // force clang-format include
                                               // order

#include <spdlog/sinks/stdout_color_sinks.h>

namespace SiFi
{

using logger = std::shared_ptr<spdlog::logger>;

inline logger createLogger(std::string name)
{
    auto alreadyExists = spdlog::get(name);
    if (alreadyExists != nullptr) { return alreadyExists; }

    auto logger = spdlog::stdout_color_st(name);
    spdlog::drop(name); // it will be released after all references are lost
    return logger;
}

} // namespace SiFi

#endif
