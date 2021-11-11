# file(GLOB cc6_cm_SRCS "*.cc")
# file(GLOB cc6_cm_HDRS "*.hh")

# ROOT_GENERATE_DICTIONARY(G__cc6_cm ${cc6_cm_HDRS} LINKDEF LinkDef.h)

# add_library(cc6_cm SHARED ${cc6_cm_SRCS} G__cc6_cm)
# target_link_libraries(cc6_cm cc6_common ${ROOT_LIBRARIES})

# set_target_properties(cc6_cm PROPERTIES
#     SOVERSION ${PROJECT_VERSION}
#     PUBLIC_HEADER "${cc6_cm_HDRS}"
# )

# install(TARGETS cc6_cm
#     LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
#     PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

CC6_GENERATE_LIBRARY(
    TARGET cc6_cm
    SOURCES
        CMSimulation.cc
        CMReconstruction.cc
        G4Reconstruction.cc
        G4SimulationAdapter.cc
        Mask.cc
    HEADERS
        CMSimulation.hh
        CMReconstruction.hh
        G4Reconstruction.hh
        G4SimulationAdapter.hh
        Mask.hh

    LIBRARIES
        ROOT::Core
        cc6_common

    PRIVATE_LIBRARIES
        cc6_tools
        spdlog::spdlog
        SiFi::CmdLineArgs
)
