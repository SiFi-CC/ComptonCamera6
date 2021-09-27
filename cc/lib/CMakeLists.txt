CC6_GENERATE_LIBRARY(
    TARGET cc6_cc
    SOURCES
        CCMLEM.cc
        CCReconstruction.cc
        CCSimulation.cc
        ComptonCone.cc
        SMElement.cc
    HEADERS
        CCMLEM.hh
        CCReconstruction.hh
        CCSimulation.hh
        ComptonCone.hh
        SMElement.hh

    LIBRARIES
        ROOT::Core

    PRIVATE_LIBRARIES
        cc6_tools
        cc6_common
        spdlog::spdlog
        SiFi::CmdLineArgs
)
