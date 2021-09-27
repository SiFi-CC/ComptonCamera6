#----------------------------------------------
# author: A.Bubak (arkadiusz@bubak.pl)
# date: 2018.01.05
# v0.01
#----------------------------------------------


#--- CMakeLists.txt for CC6 package.
cmake_minimum_required(VERSION 3.16 FATAL_ERROR)
project(ComptonCamera6
     VERSION 1.1.0
     LANGUAGES CXX
)

#--- You need to tell CMake where to find the ROOT installation.
list(APPEND CMAKE_PREFIX_PATH $ENV{ROOTSYS})
list(APPEND CMAKE_MODULE_PATH
    ${CMAKE_SOURCE_DIR}/Modules
    $ENV{ROOTSYS}
    ${CMAKE_SOURCE_DIR}/cmake-scripts
)

#---Locate the ROOT package and defines a number of variables (e.g. ROOT_INCLUDE_DIRS)
find_package(ROOT REQUIRED COMPONENTS Gpad Geom Hist MathMore Physics Tree)

find_package(CmdLineArgs 2.2.0 REQUIRED)

add_subdirectory(spdlog)

include(GNUInstallDirs)

include(helpers)

##### set paths
set(CMAKE_INSTALL_CMAKEDIR ${DEF_CMAKE_INSTALL_CMAKEDIR} CACHE PATH "Installation directory for CMake files")

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)

#---Define useful ROOT functions and macros (e.g. ROOT_GENERATE_DICTIONARY)

configure_file(macros/profile.sh.in profile.sh @ONLY)

install(FILES
    ${CMAKE_CURRENT_BINARY_DIR}/profile.sh
    macros/rootlogon.C
    DESTINATION ${CMAKE_INSTALL_DATAROOTDIR}/ComptonCamera6/)

INSTALL(DIRECTORY data/masks
    DESTINATION ${CMAKE_INSTALL_DATAROOTDIR}/ComptonCamera6/)

include(c++-standards)
include(code-coverage)
include(sanitizers)

cxx_17()
add_code_coverage()

include_directories(
    ${CMAKE_BINARY_DIR}
)

add_subdirectory(common)
add_subdirectory(tools)
add_subdirectory(cm)
add_subdirectory(cc)

option(ENABLE_TESTING "Build tests" ON)

if(ENABLE_TESTING)
    enable_testing()
    add_subdirectory(tests)
endif()

#---Locating Doxygen and generating documentation

option(BUILD_DOC "Build documentation" OFF)
option(INSTALL_DOC "Build documentation" OFF)

if (BUILD_DOC)
    find_package(Doxygen)
endif()

if(DOXYGEN_FOUND)
    set(DOXYGEN_IN ${CMAKE_CURRENT_SOURCE_DIR}/Doxyfile.in)
    set(DOXYGEN_OUT ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile)

    configure_file(${DOXYGEN_IN} ${DOXYGEN_OUT} @ONLY)
    message("Doxygen build started")

    add_custom_command(
        OUTPUT ${CMAKE_SOURCE_DIR}/html/index.html
        COMMAND ${DOXYGEN_EXECUTABLE} ${DOXYGEN_OUT}
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
        COMMENT "Generating API documentation with Doxygen"
        DEPENDS cc6_cc cc6_cm cc6_common cc6_tools
        VERBATIM)
    add_custom_target(doc ALL
        DEPENDS ${CMAKE_SOURCE_DIR}/html/index.html)

else(DOXYGEN_FOUND)
    message("Deoxygen must be installed to generate the documentation")
endif(DOXYGEN_FOUND)

if (INSTALL_DOC AND BUILD_DOC)
    INSTALL(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/html
    DESTINATION ${CMAKE_INSTALL_DATAROOTDIR}/ComptonCamera6/doc)
endif()

# messages
MESSAGE(STATUS "<<< Configuration >>>
Project:        ${CMAKE_PROJECT_NAME}

Build type      ${CMAKE_BUILD_TYPE}
Install path    ${CMAKE_INSTALL_PREFIX}

Compiler:
C               ${CMAKE_C_COMPILER}
C++             ${CMAKE_CXX_COMPILER}

Linker:
Ld              ${CMAKE_LINKER}

Compiler flags:
C               ${CMAKE_C_FLAGS}
C++             ${CMAKE_CXX_FLAGS}

Linker flags:
Executable      ${CMAKE_EXE_LINKER_FLAGS}
Module          ${CMAKE_MODULE_LINKER_FLAGS}
Shared          ${CMAKE_SHARED_LINKER_FLAGS}

Doxygen:        ${DOXYGEN_FOUND}
Build doc:      ${BUILD_DOC}
Install doc:    ${INSTALL_DOC}
\n")

foreach(p LIB BIN INCLUDE CMAKE)
    message(STATUS "CMAKE_INSTALL_${p}DIR: ${CMAKE_INSTALL_${p}DIR}")
endforeach()
