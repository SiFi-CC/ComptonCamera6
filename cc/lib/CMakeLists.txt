file(GLOB cc6_cc_SRCS "*.cc")
file(GLOB cc6_cc_HDRS "*.hh")

ROOT_GENERATE_DICTIONARY(G__cc6_cc ${cc6_cc_HDRS} LINKDEF LinkDef.h)

add_library(cc6_cc SHARED ${cc6_cc_SRCS} G__cc6_cc)
target_link_libraries(cc6_cc
    cc6_common cc6_tools
    ${ROOT_LIBRARIES})

set_target_properties(cc6_cc PROPERTIES
    SOVERSION ${PROJECT_VERSION}
    PUBLIC_HEADER "${cc6_cc_HDRS}"
)

install(TARGETS cc6_cc
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
