file(GLOB cc6_tools_SRCS "*.cc")
file(GLOB cc6_tools_HDRS "*.hh")

ROOT_GENERATE_DICTIONARY(G__cc6_tools ${cc6_tools_HDRS} LINKDEF LinkDef.h)

add_library(cc6_tools SHARED ${cc6_tools_SRCS} G__cc6_tools.cxx)
target_link_libraries(cc6_tools
    ${ROOT_LIBRARIES}
    cc6_common
)

set_target_properties(cc6_tools PROPERTIES
    SOVERSION ${PROJECT_VERSION}
    PUBLIC_HEADER "${cc6_tools_HDRS}"
)

install(TARGETS cc6_tools
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
