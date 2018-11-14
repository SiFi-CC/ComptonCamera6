add_executable(cm_reconstruct cm_reconstruct.cc)
add_executable(cm_simulation cm_simulation.cc)

target_link_libraries(cm_reconstruct cc6_cm)
target_link_libraries(cm_simulation cc6_cm)

install(TARGETS cm_reconstruct cm_simulation
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
