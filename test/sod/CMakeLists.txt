add_test(NAME sod-run COMMAND mfe1 ${CMAKE_CURRENT_SOURCE_DIR}/mfein.json)
add_test(NAME sod-cmp COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/test.py ${CMAKE_CURRENT_SOURCE_DIR}/mfegrf-ref)

set_tests_properties(sod-run PROPERTIES
    ENVIRONMENT_MODIFICATION
        LD_LIBRARY_PATH=path_list_prepend:${CMAKE_BINARY_DIR}/src
    FIXTURES_SETUP sod
)
set_tests_properties(sod-cmp PROPERTIES FIXTURES_REQUIRED sod)


