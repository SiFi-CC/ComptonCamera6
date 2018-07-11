# $Id: CMakeLists.txt Jonas Kasper $
# Purpose. Example for Simulation output Read in class $

#----------------------------------------------------------------------------
# Setup the project
cmake_minimum_required(VERSION 3.0 FATAL_ERROR)


project(Example_Project)

SET(CMAKE_INSTALL_PREFIX ${CMAKE_CURRENT_BINARY_DIR}) 
message(STATUS "Install path " ${CMAKE_INSTALL_PREFIX})

#----------Set Compiler
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} $ENV{ROOTSYS} ${CMAKE_CURRENT_SOURCE_DIR} )

message(STATUS "CMAKE_MODULE_PATH: " ${CMAKE_MODULE_PATH})
message(STATUS "Found CXXCompiler version: " ${CMAKE_CXX_COMPILER_VERSION})


#--------------Root Setup
list(APPEND CMAKE_PREFIX_PATH $ENV{ROOTSYS}) 

message(STATUS "CMAKE_PREFIX_PATH: " ${CMAKE_PREFIX_PATH})
find_package( ROOT REQUIRED COMPONENTS Gui Core) 
message( STATUS "Found ROOT version " ${ROOT_VERSION} ) 

include(${ROOT_USE_FILE})

include(makeRootDict)


include_directories(${ROOT_INCLUDE_DIRS} )

# Setup include directory for this project

include_directories(${PROJECT_SOURCE_DIR}/include)

#---------------------------------------------------------------------------
#create extra libraries

MAKE_ROOT_DICT(DICTNAME StructLib LINKDEF include/StructLinkDef.h SOURCES src/DR_GenerallStructs.cc INCLUDES include/Structinclude.hh include/DR_GenerallStructs.hh LIBRARIES ${ROOT_LIBRARIES})

#----------------------------------------------------------------------------
# Locate sources and headers for this project
# NB: headers are included so they will show up in IDEs


file(GLOB sources ${PROJECT_SOURCE_DIR}/src/*.cc)
file(GLOB headers ${PROJECT_SOURCE_DIR}/include/*.hh)

#----------------------------------------------------------------------------
# Add the executable, and link it to the libraries


SET(CMAKE_INSTALL_RPATH ${CMAKE_INSTALL_PREFIX}/lib )


add_executable(ReadIn ReadIn.cc ${sources} ${headers})
target_link_libraries(ReadIn StructLib ${ROOT_LIBRARIES})

#----------------------------------------------------------------------------

add_custom_target(Example_Project DEPENDS ReadIn)

#----------------------------------------------------------------------------
# Install the executable to 'bin' directory under CMAKE_INSTALL_PREFIX

set_property(TARGET ReadIn PROPERTY INSTALL_RPATH_USE_LINK_PATH TRUE)
install(TARGETS ReadIn RUNTIME DESTINATION ${CMAKE_INSTALL_PREFIX}/bin)

