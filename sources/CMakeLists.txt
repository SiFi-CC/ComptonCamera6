#----------------------------------------------
# author: A.Bubak (arkadiusz@bubak.pl)
# date: 2018.01.05
# v0.01
#----------------------------------------------


#--- CMakeLists.txt for CC6 package.
cmake_minimum_required(VERSION 3.0 FATAL_ERROR)
project(ComptonCamera6)


#--- You need to tell CMake where to find the ROOT installation.
list(APPEND CMAKE_PREFIX_PATH $ENV{ROOTSYS})


#---Locate the ROOT package and defines a number of variables (e.g. ROOT_INCLUDE_DIRS)
find_package(ROOT REQUIRED COMPONENTS RIO)


#---Set global variables e.g. where sources and includes are
file(GLOB sources ${PROJECT_SOURCE_DIR}/src/*.cc)
file(GLOB headers ${PROJECT_SOURCE_DIR}/include/*.hh)


#---Define useful ROOT functions and macros (e.g. ROOT_GENERATE_DICTIONARY)
include(${ROOT_USE_FILE})
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include)
ROOT_GENERATE_DICTIONARY(G__ComptonCamera6 ${headers} LINKDEF LinkDef.hh)


#---Create a shared library with generated dictionary
add_library(ComptonCamera6 SHARED camera.cc ${sources} G__ComptonCamera6.cxx)
target_link_libraries(ComptonCamera6 ${ROOT_LIBRARIES})


#---Create  a main program using the library
add_executable(camera camera.cc)
add_executable(mask mask.cc)
add_executable(reconstruction reconstruction.cc)
add_executable(cmreco cmreco.cc)

target_link_libraries(camera ComptonCamera6)
target_link_libraries(mask ComptonCamera6)
target_link_libraries(reconstruction ComptonCamera6)
target_link_libraries(cmreco ComptonCamera6)
