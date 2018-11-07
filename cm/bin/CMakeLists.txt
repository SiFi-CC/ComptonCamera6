add_executable(cmreco cmreco.cc)

target_link_libraries(cmreco cc6_cm)

install(TARGETS cmreco
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
