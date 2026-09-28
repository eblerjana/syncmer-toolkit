add_library(syncmer-lib STATIC
	commands.cpp
	syncmerstats.cpp
)

target_include_directories(syncmer-lib PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(syncmer-lib PUBLIC syngcore)
target_compile_options(syncmer-lib PRIVATE -Wno-write-strings)

add_executable(syncmer-toolkit
	main.cpp
)

target_link_libraries(syncmer-toolkit PRIVATE syncmer-lib)
target_link_options(syncmer-toolkit PRIVATE -static-libgcc -static-libstdc++)
target_compile_options(syncmer-toolkit PRIVATE -Wno-write-strings)
