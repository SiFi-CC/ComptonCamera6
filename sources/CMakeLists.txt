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
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} $ENV{ROOTSYS} ${CMAKE_CURRENT_SOURCE_DIR} )

#---Locate the ROOT package and defines a number of variables (e.g. ROOT_INCLUDE_DIRS)
find_package(ROOT REQUIRED COMPONENTS RIO Geom)


#---Set global variables e.g. where sources and includes are
file(GLOB sources ${PROJECT_SOURCE_DIR}/src/*.cc)
file(GLOB headers ${PROJECT_SOURCE_DIR}/include/*.hh)


#---Define useful ROOT functions and macros (e.g. ROOT_GENERATE_DICTIONARY)
include(makeRootDict)
include(${ROOT_USE_FILE})
include_directories(${ROOT_INCLUDE_DIRS})
include_directories(${PROJECT_SOURCE_DIR}/include)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include)

ROOT_GENERATE_DICTIONARY(G__ComptonCamera6 ${headers} LINKDEF LinkDef.hh)

MAKE_ROOT_DICT(DICTNAME StructLib LINKDEF StructLinkDef.hh SOURCES src/DR_GenerallStructs.cc INCLUDES include/Structinclude.hh include/DR_GenerallStructs.hh LIBRARIES ${ROOT_LIBRARIES})

#---Create a shared library with generated dictionary
add_library(ComptonCamera6 SHARED camera.cc ${sources} G__ComptonCamera6.cxx)
target_link_libraries(ComptonCamera6 ${ROOT_LIBRARIES} StructLib)


#---Create  a main program using the library
add_executable(camera camera.cc)
add_executable(mask mask.cc)
add_executable(reconstruction reconstruction.cc)
add_executable(cmreco cmreco.cc)
add_executable(mlem mlem.cc)
add_executable(tests tests.cc)

target_link_libraries(camera ComptonCamera6)
target_link_libraries(mask ComptonCamera6)
target_link_libraries(reconstruction ComptonCamera6)
target_link_libraries(cmreco ComptonCamera6)
target_link_libraries(mlem ComptonCamera6)
target_link_libraries(tests ComptonCamera6)


#---Locating Doxygen and generating documentation 
find_package(Doxygen)

if(DOXYGEN_FOUND)
    set(DOXYGEN_IN ${CMAKE_CURRENT_SOURCE_DIR}/../sources/Doxyfile.in)
    set(DOXYGEN_OUT ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile)

    configure_file(${DOXYGEN_IN} ${DOXYGEN_OUT} @ONLY)
    message("Doxygen build started")
    
    add_custom_target(doc_doxygen ALL
        COMMAND ${DOXYGEN_EXECUTABLE} ${DOXYGEN_OUT}
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
        COMMENT "Generating API documentation with Doxygen"
        VERBATIM)
else(DOXYGEN_FOUND)
    message("Deoxygen must be installed to generate the documentation")
endif(DOXYGEN_FOUND)