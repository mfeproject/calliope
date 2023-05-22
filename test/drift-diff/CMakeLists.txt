add_test(NAME drift-diff-run COMMAND mfe1 ${CMAKE_CURRENT_SOURCE_DIR}/mfein.json)
add_test(NAME drift-diff-cmp COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/test.py ${CMAKE_CURRENT_SOURCE_DIR}/mfegrf-ref)

set_tests_properties(drift-diff-run PROPERTIES
    ENVIRONMENT_MODIFICATION
        LD_LIBRARY_PATH=path_list_prepend:${CMAKE_BINARY_DIR}/src
    FIXTURES_SETUP drift-diff
)
set_tests_properties(drift-diff-cmp PROPERTIES FIXTURES_REQUIRED drift-diff)


