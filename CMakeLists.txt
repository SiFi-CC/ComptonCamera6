#----------------------------------------------
# author: A.Bubak (arkadiusz@bubak.pl)
# date: 2018.01.05
# v0.01
#----------------------------------------------


#--- CMakeLists.txt for CC6 package.
cmake_minimum_required(VERSION 3.0 FATAL_ERROR)
project(ComptonCamera6 CXX)
string(TOUPPER "${PROJECT_NAME}" CAPS_PROJECT_NAME)

set (VERSION_MAJOR 0)
set (VERSION_MINOR 0)
set (VERSION_PATCH 1)
set(PROJECT_VERSION ${VERSION_MAJOR}.${VERSION_MINOR}.${VERSION_PATCH})

#--- You need to tell CMake where to find the ROOT installation.
list(APPEND CMAKE_PREFIX_PATH $ENV{ROOTSYS})
list(APPEND CMAKE_MODULE_PATH
    ${CMAKE_SOURCE_DIR}/Modules
    $ENV{ROOTSYS}
    ${CMAKE_SOURCE_DIR}/cmake-scripts
)

#---Locate the ROOT package and defines a number of variables (e.g. ROOT_INCLUDE_DIRS)
find_package(ROOT REQUIRED COMPONENTS RIO Geom MathMore)
set(ROOTSYS $ENV{ROOTSYS})

find_package(CmdLineArgs REQUIRED)

add_subdirectory(spdlog)
include_directories(spdlog/include)

include(GNUInstallDirs)

##### set paths
if(WIN32 AND NOT CYGWIN)
    set(DEF_CMAKE_INSTALL_CMAKEDIR cmake)
else()
    set(DEF_CMAKE_INSTALL_CMAKEDIR ${CMAKE_INSTALL_LIBDIR}/cmake/${CMAKE_PROJECT_NAME})
endif()
set(CMAKE_INSTALL_CMAKEDIR ${DEF_CMAKE_INSTALL_CMAKEDIR} CACHE PATH "Installation directory for CMake files")

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)

#---Define useful ROOT functions and macros (e.g. ROOT_GENERATE_DICTIONARY)

include(${ROOT_USE_FILE})
include_directories(${ROOT_INCLUDE_DIRS})

configure_file(macros/profile.sh.in profile.sh @ONLY)

install(FILES
    ${CMAKE_CURRENT_BINARY_DIR}/profile.sh
    macros/rootlogon.C
    DESTINATION ${CMAKE_INSTALL_DATAROOTDIR}/ComptonCamera6/)

INSTALL(DIRECTORY data/masks
    DESTINATION ${CMAKE_INSTALL_DATAROOTDIR}/ComptonCamera6/)

include_directories(common/lib)
include_directories(tools/lib)

link_directories(${ROOT_LIBRARY_DIRS})

include(c++-standards)
include(code-coverage)
include(sanitizers)

cxx_11()
add_code_coverage()

add_subdirectory(common)
add_subdirectory(tools)
add_subdirectory(cm)
add_subdirectory(cc)

enable_testing()
add_subdirectory(tests)

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

# unset  (ELOG CACHE)
if (ELOG)
# if you are building in-source, this is the same as CMAKE_SOURCE_DIR, otherwise
# this is the top level directory of your build tree
MESSAGE( STATUS "CMAKE_BINARY_DIR:         " ${CMAKE_BINARY_DIR} )
# if you are building in-source, this is the same as CMAKE_CURRENT_SOURCE_DIR, otherwise this
# is the directory where the compiled or generated files from the current CMakeLists.txt will go to
MESSAGE( STATUS "CMAKE_CURRENT_BINARY_DIR: " ${CMAKE_CURRENT_BINARY_DIR} )
# this is the directory, from which cmake was started, i.e. the top level source directory
MESSAGE( STATUS "CMAKE_SOURCE_DIR:         " ${CMAKE_SOURCE_DIR} )
# this is the directory where the currently processed CMakeLists.txt is located in
MESSAGE( STATUS "CMAKE_CURRENT_SOURCE_DIR: " ${CMAKE_CURRENT_SOURCE_DIR} )
# contains the full path to the top level directory of your build tree
MESSAGE( STATUS "PROJECT_BINARY_DIR: " ${PROJECT_BINARY_DIR} )
# contains the full path to the root of your project source directory,
# i.e. to the nearest directory where CMakeLists.txt contains the PROJECT() command
MESSAGE( STATUS "PROJECT_SOURCE_DIR: " ${PROJECT_SOURCE_DIR} )
# set this variable to specify a common place where CMake should put all executable files
# (instead of CMAKE_CURRENT_BINARY_DIR)
MESSAGE( STATUS "EXECUTABLE_OUTPUT_PATH: " ${EXECUTABLE_OUTPUT_PATH} )
# set this variable to specify a common place where CMake should put all libraries
# (instead of CMAKE_CURRENT_BINARY_DIR)
MESSAGE( STATUS "LIBRARY_OUTPUT_PATH:     " ${LIBRARY_OUTPUT_PATH} )
# tell CMake to search first in directories listed in CMAKE_MODULE_PATH
# when you use FIND_PACKAGE() or INCLUDE()
MESSAGE( STATUS "CMAKE_MODULE_PATH: " ${CMAKE_MODULE_PATH} )
# this is the complete path of the cmake which runs currently (e.g. /usr/local/bin/cmake)
MESSAGE( STATUS "CMAKE_COMMAND: " ${CMAKE_COMMAND} )
# this is the CMake installation directory
MESSAGE( STATUS "CMAKE_ROOT: " ${CMAKE_ROOT} )
# this is the filename including the complete path of the file where this variable is used.
MESSAGE( STATUS "CMAKE_CURRENT_LIST_FILE: " ${CMAKE_CURRENT_LIST_FILE} )
# this is linenumber where the variable is used
MESSAGE( STATUS "CMAKE_CURRENT_LIST_LINE: " ${CMAKE_CURRENT_LIST_LINE} )
# this is used when searching for include files e.g. using the FIND_PATH() command.
MESSAGE( STATUS "CMAKE_INCLUDE_PATH: " ${CMAKE_INCLUDE_PATH} )
# this is used when searching for libraries e.g. using the FIND_LIBRARY() command.
MESSAGE( STATUS "CMAKE_LIBRARY_PATH: " ${CMAKE_LIBRARY_PATH} )
# the complete system name, e.g. "Linux-2.4.22", "FreeBSD-5.4-RELEASE" or "Windows 5.1"
MESSAGE( STATUS "CMAKE_SYSTEM: " ${CMAKE_SYSTEM} )
# the short system name, e.g. "Linux", "FreeBSD" or "Windows"
MESSAGE( STATUS "CMAKE_SYSTEM_NAME: " ${CMAKE_SYSTEM_NAME} )
# only the version part of CMAKE_SYSTEM
MESSAGE( STATUS "CMAKE_SYSTEM_VERSION: " ${CMAKE_SYSTEM_VERSION} )
# the processor name (e.g. "Intel(R) Pentium(R) M processor 2.00GHz")
MESSAGE( STATUS "CMAKE_SYSTEM_PROCESSOR: " ${CMAKE_SYSTEM_PROCESSOR} )
# is TRUE on all UNIX-like OS's, including Apple OS X and CygWin
MESSAGE( STATUS "UNIX: " ${UNIX} )
# is TRUE on Windows, including CygWin
MESSAGE( STATUS "WIN32: " ${WIN32} )
# is TRUE on Apple OS X
MESSAGE( STATUS "APPLE: " ${APPLE} )
# is TRUE when using the MinGW compiler in Windows
MESSAGE( STATUS "MINGW: " ${MINGW} )
# is TRUE on Windows when using the CygWin version of cmake
MESSAGE( STATUS "CYGWIN: " ${CYGWIN} )
# is TRUE on Windows when using a Borland compiler
MESSAGE( STATUS "BORLAND: " ${BORLAND} )
# Microsoft compiler
MESSAGE( STATUS "MSVC: " ${MSVC} )
MESSAGE( STATUS "MSVC_IDE: " ${MSVC_IDE} )
MESSAGE( STATUS "MSVC60: " ${MSVC60} )
MESSAGE( STATUS "MSVC70: " ${MSVC70} )
MESSAGE( STATUS "MSVC71: " ${MSVC71} )
MESSAGE( STATUS "MSVC80: " ${MSVC80} )
MESSAGE( STATUS "CMAKE_COMPILER_2005: " ${CMAKE_COMPILER_2005} )
# set this to true if you don't want to rebuild the object files if the rules have changed,
# but not the actual source files or headers (e.g. if you changed the some compiler switches)
MESSAGE( STATUS "CMAKE_SKIP_RULE_DEPENDENCY: " ${CMAKE_SKIP_RULE_DEPENDENCY} )
# since CMake 2.1 the install rule depends on all, i.e. everything will be built before installing.
# If you don't like this, set this one to true.
MESSAGE( STATUS "CMAKE_SKIP_INSTALL_ALL_DEPENDENCY: " ${CMAKE_SKIP_INSTALL_ALL_DEPENDENCY} )
# If set, runtime paths are not added when using shared libraries. Default it is set to OFF
MESSAGE( STATUS "CMAKE_SKIP_RPATH: " ${CMAKE_SKIP_RPATH} )
# set this to true if you are using makefiles and want to see the full compile and link
# commands instead of only the shortened ones
MESSAGE( STATUS "CMAKE_VERBOSE_MAKEFILE: " ${CMAKE_VERBOSE_MAKEFILE} )
# this will cause CMake to not put in the rules that re-run CMake. This might be useful if
# you want to use the generated build files on another machine.
MESSAGE( STATUS "CMAKE_SUPPRESS_REGENERATION: " ${CMAKE_SUPPRESS_REGENERATION} )
# A simple way to get switches to the compiler is to use ADD_DEFINITIONS().
# But there are also two variables exactly for this purpose:
# the compiler flags for compiling C sources
MESSAGE( STATUS "CMAKE_C_FLAGS: " ${CMAKE_C_FLAGS} )
# the compiler flags for compiling C++ sources
MESSAGE( STATUS "CMAKE_CXX_FLAGS: " ${CMAKE_CXX_FLAGS} )
# Choose the type of build.  Example: SET(CMAKE_BUILD_TYPE Debug)
MESSAGE( STATUS "CMAKE_BUILD_TYPE: " ${CMAKE_BUILD_TYPE} )
# if this is set to ON, then all libraries are built as shared libraries by default.
MESSAGE( STATUS "BUILD_SHARED_LIBS: " ${BUILD_SHARED_LIBS} )
# the compiler used for C files
MESSAGE( STATUS "CMAKE_C_COMPILER: " ${CMAKE_C_COMPILER} )
# the compiler used for C++ files
MESSAGE( STATUS "CMAKE_CXX_COMPILER: " ${CMAKE_CXX_COMPILER} )
# if the compiler is a variant of gcc, this should be set to 1
MESSAGE( STATUS "CMAKE_COMPILER_IS_GNUCC: " ${CMAKE_COMPILER_IS_GNUCC} )
# if the compiler is a variant of g++, this should be set to 1
MESSAGE( STATUS "CMAKE_COMPILER_IS_GNUCXX : " ${CMAKE_COMPILER_IS_GNUCXX} )
# the tools for creating libraries
MESSAGE( STATUS "CMAKE_AR: " ${CMAKE_AR} )
MESSAGE( STATUS "CMAKE_RANLIB: " ${CMAKE_RANLIB} )

endif (ELOG)
