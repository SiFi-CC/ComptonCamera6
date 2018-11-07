add_executable(camera camera.cc)
add_executable(mlem mlem.cc)
add_executable(reconstruction reconstruction.cc)

target_link_libraries(camera cc6_cc)
target_link_libraries(mlem cc6_cc cc6_tools)
target_link_libraries(reconstruction cc6_cc cc6_common cc6_tools)

install(TARGETS reconstruction mlem camera
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
