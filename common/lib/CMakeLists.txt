file(GLOB cc6_common_SRCS "*.cc" "**/*.cc")
file(GLOB cc6_common_HDRS "*.hh" "**/*.hh")

ROOT_GENERATE_DICTIONARY(G__cc6_common ${cc6_common_HDRS} LINKDEF LinkDef.h)

add_library(cc6_common SHARED ${cc6_common_SRCS} G__cc6_common.cxx)
target_link_libraries(cc6_common cc6_tools
    ${ROOT_LIBRARIES})

set_target_properties(cc6_common PROPERTIES
    SOVERSION ${PROJECT_VERSION}
    PUBLIC_HEADER "${cc6_common_HDRS}"
)

install(TARGETS cc6_common
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
