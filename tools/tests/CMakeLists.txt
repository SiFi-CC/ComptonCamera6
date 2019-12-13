set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR})
find_package(cppunit)

file(GLOB UnitTests_SRCS RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} "*_tests.cpp" )

link_directories(${CMAKE_CURRENT_BINARY_DIR}/lib)

add_executable(UnitTesterTools test_runner.cpp ${UnitTests_SRCS})
target_link_libraries(UnitTesterTools CmdLineArgs ${CPPUNIT_LIBRARIES} cc6_common cc6_tools)

foreach(test ${UnitTests_SRCS})
    file(READ ${test} test_file)
    string(REGEX MATCH "CPPUNIT_TEST_SUITE\(.([A-Za-z0-9_ ]+).\)" _ ${test_file})
    string(STRIP ${CMAKE_MATCH_2} suite_name)
    message(STATUS "Add test ${suite_name}")
    add_test(NAME ${suite_name} COMMAND UnitTesterTools ${suite_name})
endforeach(test)
