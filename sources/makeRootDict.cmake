include(CMakeParseArguments)
function(MAKE_ROOT_DICT)
 	set(options NOINSTALL)
 	set(oneValueArgs DICTNAME LINKDEF)
 	set(multiValueArgs SOURCES INCLUDES LIBRARIES)
 	cmake_parse_arguments(MAKE_ROOT_DICT "${options}" "${oneValueArgs}" "${multiValueArgs}" ${ARGN} ) 
	if(${ROOT_VERSION_MAJOR} EQUAL 6 ) 
		root_generate_dictionary(G__${MAKE_ROOT_DICT_DICTNAME} ${MAKE_ROOT_DICT_INCLUDES} MODULE ${MAKE_ROOT_DICT_DICTNAME} LINKDEF ${MAKE_ROOT_DICT_LINKDEF})
 	else()
 		root_generate_dictionary(G__${MAKE_ROOT_DICT_DICTNAME} ${MAKE_ROOT_DICT_INCLUDES} LINKDEF ${MAKE_ROOT_DICT_LINKDEF} OPTIONS -p)
 	endif()
 	add_library(${MAKE_ROOT_DICT_DICTNAME} SHARED ${MAKE_ROOT_DICT_SOURCES} G__${MAKE_ROOT_DICT_DICTNAME}.cxx)
 	target_link_libraries(${MAKE_ROOT_DICT_DICTNAME} ${MAKE_ROOT_DICT_LIBRARIES})
 	target_include_directories ( ${MAKE_ROOT_DICT_DICTNAME} PUBLIC ${CMAKE_CURRENT_SOURCE_DIR} )
 	if(MAKE_ROOT_DICT_NOINSTALL)
 	else()
 		install(TARGETS ${MAKE_ROOT_DICT_DICTNAME} DESTINATION lib)
			if(${ROOT_VERSION_MAJOR} EQUAL 6 )
 				install(FILES
 					${CMAKE_CURRENT_BINARY_DIR}/lib${MAKE_ROOT_DICT_DICTNAME}_rdict.pcm
 					${CMAKE_CURRENT_BINARY_DIR}/lib${MAKE_ROOT_DICT_DICTNAME}.rootmap
 					DESTINATION lib)
 			endif()
 	endif()
endfunction(MAKE_ROOT_DICT)


