set(SYNG_DIR ${CMAKE_SOURCE_DIR}/external/syng)

option(SYNG_BAMIO "enable SAM/BAM/CRAM support (requires htslib)" OFF)

set(SYNG_CORE_SOURCES
	${SYNG_DIR}/hash.c
	${SYNG_DIR}/dict.c
	${SYNG_DIR}/array.c
	${SYNG_DIR}/utils.c
	${SYNG_DIR}/ONElib.c
	${SYNG_DIR}/seqio.c
	${SYNG_DIR}/seqhash.c
	${SYNG_DIR}/kmerhash.c
	${SYNG_DIR}/rskip.c
	${SYNG_DIR}/syncmerset.c
	${SYNG_DIR}/syngbwt3.c
)

add_library(syngcore STATIC ${SYNG_CORE_SOURCES})
target_include_directories(syngcore PUBLIC ${SYNG_DIR})

target_compile_definitions(syngcore PRIVATE
	ONEIO
	$<$<BOOL:${SYNG_BAMIO}>:BAMIO>
)

find_package(ZLIB REQUIRED)
find_package(Threads REQUIRED)

target_link_libraries(syngcore PUBLIC ZLIB::ZLIB Threads::Threads m)

if(SYNG_BAMIO)
	set(HTSLIB_DIR "${CMAKE_SOURCE_DIR}/external/htslib" CACHE PATH "Path to htslib")
	target_include_directories(syngcore PUBLIC ${HTSLIB_DIR})
	find_library(HTS_LIB hts PATHS ${HTSLIB_DIR} REQUIRED)
	target_link_libraries(syngcore PUBLIC ${HTS_LIB} bz2 lzma curl)
endif()
