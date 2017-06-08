
SET(FCLIB_TEST_PROGRAMS
	fctst_cpp
	fclib_matrix_test
)

foreach(TEST_SEL ${FCLIB_TEST_PROGRAMS})
	add_executable (${TEST_SEL}
					${TEST_SEL}
					"${FCLIB_CPP_INCLUDEDIR}/${FCLIB_CPP_HEADERS}")
					
	target_include_directories (${TEST_SEL}
					PUBLIC ${FCLIB_CPP_INCLUDEDIR})
					
	target_link_libraries(${TEST_SEL} PUBLIC ${FCLIB_CPP_LIBS})

endforeach(TEST_SEL)
