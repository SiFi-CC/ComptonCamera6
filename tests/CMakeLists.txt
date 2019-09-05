set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR})
find_package(cppunit)

file(GLOB UnitTests_SRCS RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} "*_tests.cpp" )

link_directories(${CMAKE_CURRENT_BINARY_DIR}/lib)

add_executable(UnitTester test_runner.cpp ${UnitTests_SRCS})
target_link_libraries(UnitTester CmdLineArgs ${CPPUNIT_LIBRARIES})

foreach(test ${UnitTests_SRCS})
    file(READ ${test} test_file)
    string(REGEX MATCH "CPPUNIT_TEST_SUITE\(.([A-Za-z0-9 ]+).\)" _ ${test_file})
    string(STRIP ${CMAKE_MATCH_2} suite_name)
    add_test(NAME ${suite_name} COMMAND UnitTester ${suite_name})
endforeach(test)
