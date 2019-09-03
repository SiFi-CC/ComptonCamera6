add_executable(cm_reconstruct cm_reconstruct.cc)
add_executable(cm_simulation cm_simulation.cc)
add_executable(cm_g4reconstruct cm_g4reconstruct.cc)

target_link_libraries(cm_reconstruct cc6_cm CmdLineArgs)
target_link_libraries(cm_simulation cc6_cm CmdLineArgs)
target_link_libraries(cm_g4reconstruct cc6_cm CmdLineArgs)

install(TARGETS cm_reconstruct cm_simulation cm_g4reconstruct
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
