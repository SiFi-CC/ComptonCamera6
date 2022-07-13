FetchContent_Declare(spdlog
    GIT_REPOSITORY https://github.com/gabime/spdlog
    GIT_TAG        v${REQUIRED_SPDLOG_VERSION}
)

if(SIFI_BUILTIN_SPDLOG STREQUAL "AUTO")
    find_package(spdlog ${REQUIRED_SPDLOG_VERSION} QUIET)
    if (NOT spdlog_FOUND)
        SET(USE_BUILTIN_SPDLOG TRUE)
    endif()
elseif(SIFI_BUILTIN_SPDLOG)     # a true value (such as ON) was used
    SET(USE_BUILTIN_SPDLOG TRUE)
else()                  # a false value (such as OFF) was used
    find_package(spdlog ${REQUIRED_SPDLOG_VERSION} REQUIRED)
endif()

if (USE_BUILTIN_SPDLOG)
    #FetchContent_MakeAvailable(spdlog)
    FetchContent_GetProperties(spdlog)
    if(NOT spdlog_POPULATED)
        FetchContent_Populate(spdlog)
        add_subdirectory(${spdlog_SOURCE_DIR} ${spdlog_BINARY_DIR} EXCLUDE_FROM_ALL)
    endif()
    set_property(TARGET spdlog PROPERTY POSITION_INDEPENDENT_CODE ON)
else()
    message(STATUS "Uses system-provided spdlog")
endif()
