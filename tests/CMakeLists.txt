add_executable(unit_tests test_main.cpp test_syncmerstats.cpp test_commands.cpp)
target_include_directories(unit_tests PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/third_party)
target_link_libraries(unit_tests PRIVATE syncmer-lib)
target_compile_definitions(unit_tests PRIVATE
	TEST_DATA_DIR="${CMAKE_CURRENT_SOURCE_DIR}/data/"
	TEST_OUTPUT_DIR="${CMAKE_CURRENT_BINARY_DIR}/"
)
add_test(NAME unit_tests COMMAND unit_tests)
