file(GLOB cc6_cm_SRCS "*.cc")
file(GLOB cc6_cm_HDRS "*.hh")

ROOT_GENERATE_DICTIONARY(G__cc6_cm ${cc6_cm_HDRS} LINKDEF LinkDef.h)

add_library(cc6_cm SHARED ${cc6_cm_SRCS} G__cc6_cm)
target_link_libraries(cc6_cm cc6_common cc6_tools ${ROOT_LIBRARIES})

set_target_properties(cc6_cm PROPERTIES
    SOVERSION ${PROJECT_VERSION}
    PUBLIC_HEADER "${cc6_cm_HDRS}"
)

install(TARGETS cc6_cm
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
