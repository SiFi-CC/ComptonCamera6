find_package(GTest)

configure_file(tests_config.h.in tests_config.h)

#file(GLOB tests_SRCS RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} "tests_*.cc")
set(tests_SRCS
    tests_PhysicsBase.cc
)

add_executable(
    gtests
    ${tests_SRCS}
)

if (TARGET GTest::gtest_main)
    # new cmake provided
    set(GTEST_TRG GTest::gtest_main)
else()
    # old cmake and pc files
    set(GTEST_TRG gtest gtest_main)
endif()

target_link_libraries(
    gtests
    PRIVATE
        ${GTEST_TRG}
        #cc6_tools
        #cc6_common
        #cc6_cc
        #cc6_cm
)
#target_code_coverage(gtests ALL)

gtest_discover_tests(gtests)
