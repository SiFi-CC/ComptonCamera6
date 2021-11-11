CC6_GENERATE_LIBRARY(
    TARGET cc6_tools
    SOURCES
        DataStructConvert.cc
        PhysicsBase.cc
        Smoothing.cc

    HEADERS
        CLog.hh
        Coordinates.hh
        DataStructConvert.hh
        PhysicsBase.hh
        Smoothing.hh

    LIBRARIES
        spdlog::spdlog
        ROOT::Core

    PRIVATE_LIBRARIES
        ROOT::Matrix
        ROOT::Physics
        ROOT::Hist
)
