cmake_minimum_required(VERSION 3.16)
project(melon3d VERSION 1.2.1 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# When consumed via add_subdirectory / FetchContent, build just the library
if(NOT DEFINED PROJECT_IS_TOP_LEVEL)
	if(CMAKE_CURRENT_SOURCE_DIR STREQUAL CMAKE_SOURCE_DIR)
		set(PROJECT_IS_TOP_LEVEL TRUE)
	else()
		set(PROJECT_IS_TOP_LEVEL FALSE)
	endif()
endif()

if(NOT CMAKE_BUILD_TYPE)
	set(CMAKE_BUILD_TYPE Release)
endif()

add_library(melon3d STATIC
	src/collide.cpp
	src/grid.cpp
	src/xpbd.cpp
	src/world.cpp
)
add_library(melon3d::melon3d ALIAS melon3d)
target_include_directories(melon3d PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)

# SIMD flags only where the target actually is x86; Apple Silicon and other
# ARM builds still autovectorize the branchless lane loops with NEON
set(MELON3D_X86 FALSE)
if(CMAKE_SYSTEM_PROCESSOR MATCHES "^(x86_64|AMD64|amd64|i[3-6]86)$")
	set(MELON3D_X86 TRUE)
endif()

if(MSVC)
	if(MELON3D_X86)
		target_compile_options(melon3d PRIVATE $<$<CONFIG:Release>:/arch:AVX2>)
	endif()
else()
	target_compile_options(melon3d PRIVATE -Wall -Wextra)
	target_compile_options(melon3d PRIVATE $<$<CONFIG:Release>:-O3>)
	if(MELON3D_X86)
		target_compile_options(melon3d PRIVATE $<$<CONFIG:Release>:-mavx2 -mfma>)
	endif()
endif()

find_package(Threads REQUIRED)
target_link_libraries(melon3d PUBLIC Threads::Threads)

option(MELON3D_BUILD_TESTS "Build Melon3D tests" ${PROJECT_IS_TOP_LEVEL})
option(MELON3D_BUILD_BENCHMARKS "Build Melon3D benchmarks" ${PROJECT_IS_TOP_LEVEL})

if(MELON3D_BUILD_TESTS)
	enable_testing()
	add_executable(melon3d_tests test/test_main.cpp)
	target_link_libraries(melon3d_tests PRIVATE melon3d)
	add_test(NAME melon3d_tests COMMAND melon3d_tests)
endif()

if(MELON3D_BUILD_BENCHMARKS)
	add_executable(melon3d_bench bench/bench.cpp)
	target_link_libraries(melon3d_bench PRIVATE melon3d)
endif()

include(GNUInstallDirs)
install(TARGETS melon3d ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR})
install(FILES include/melon3d.h DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
