add_executable(syncmer-toolkit main.cpp)
target_link_libraries(syncmer-toolkit PRIVATE syngcore)
target_compile_options(syncmer-toolkit PRIVATE -Wno-write-strings)
