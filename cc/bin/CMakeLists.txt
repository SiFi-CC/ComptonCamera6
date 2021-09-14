add_executable(cc_simulation cc_simulation.cc)
add_executable(mlem mlem.cc)
add_executable(cc_reconstruct cc_reconstruct.cc)
add_executable(smatrix smatrix.cc)

target_link_libraries(cc_simulation cc6_cc CmdLineArgs)
target_link_libraries(mlem cc6_cc cc6_tools CmdLineArgs)
target_link_libraries(cc_reconstruct cc6_cc cc6_common cc6_tools CmdLineArgs)
target_link_libraries(smatrix cc6_cc cc6_tools CmdLineArgs)

install(TARGETS cc_simulation cc_reconstruct mlem smatrix
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})

