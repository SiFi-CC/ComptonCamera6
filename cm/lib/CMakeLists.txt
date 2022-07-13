CC6_GENERATE_LIBRARY(
    TARGET cc6_cm
    SOURCES
        CMReconstruction.cc
        CMSimulation.cc
        G4Reconstruction.cc
        G4SimulationAdapter.cc
        Mask.cc

    HEADERS
        CMReconstruction.hh
        CMSimulation.hh
        G4Reconstruction.hh
        G4SimulationAdapter.hh
        Mask.hh

    LIBRARIES
        cc6_common

    PRIVATE_LIBRARIES
        cc6_tools
        ROOT::Core
        spdlog::spdlog
        SiFi::CmdLineArgs
)
