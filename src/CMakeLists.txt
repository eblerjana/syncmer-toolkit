add_executable(syncmer-toolkit
	main.cpp
	syncmerstats.cpp
)
target_link_libraries(syncmer-toolkit PRIVATE syngcore)
target_link_options(syncmer-toolkit PRIVATE -static-libgcc -static-libstdc++)
target_compile_options(syncmer-toolkit PRIVATE -Wno-write-strings)
