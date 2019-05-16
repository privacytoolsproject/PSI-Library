cmake_minimum_required (VERSION 3.12.4)
project (dpfc CXX)
set (dpfc_VERSION_MAJOR 1)
set (dpfc_VERSION_MINOR 0)

# add compiler flags
if (NOT CMAKE_COMPILER_IS_GNUCXX)
  set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} --std=c++17 -g -Wall -Wextra -Wfloat-equal -Wundef -Wcas\
t-align -Wwrite-strings -Wlogical-op -Wmissing-declarations -Wredundant-decls -Wshadow -Woverloaded\
-virtual -pedantic")
endif ()

# set output directories
set (CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set (CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

set (SRCS utils count hist_sparse hist_compact)
set (TSRCS utils count hist_simple hist_sparse hist_compact)

# find required dependencies
find_package(GTest REQUIRED)

# TODO write a find_package for gmp, ntl and sodium
set (GMPXX_LIBRARIES gmpxx gmp)
set (NTL_LIBRARIES ntl)
set (SODIUM_LIBRARIES sodium)

set (LIBRARIES ntl gmpxx gmp sodium ${GTEST_LIBRARIES})
include_directories (/c/msys64/mingw64/include ${CMAKE_SOURCE_DIR}/include)

# create library
set (SRCS_LIB ${SRCS})
list (TRANSFORM SRCS_LIB PREPEND ${CMAKE_SOURCE_DIR}/src/)

add_library (dpfc ${SRCS_LIB})

link_directories (CMAKE_ARCHIVE_OUTPUT_DIRECTORY/dpfc)

# create tests
enable_testing ()

list (TRANSFORM TSRCS PREPEND t_)

foreach (SRC IN LISTS TSRCS)
	add_executable (${SRC} ${CMAKE_SOURCE_DIR}/test/${SRC})
	target_include_directories (${SRC} PUBLIC ${GTEST_INCLUDE_DIRS})
	target_link_libraries (${SRC} dpfc ${LIBRARIES})
	gtest_discover_tests (${SRC})
endforeach (SRC)

# create demos

set (DEMOS d_benchmark d_string)

foreach (DEMO IN LISTS DEMOS)
	add_executable (${DEMO} ${CMAKE_SOURCE_DIR}/demos/${DEMO})
	target_link_libraries (${DEMO} dpfc ${LIBRARIES})
	target_include_directories (${DEMO} PUBLIC ${CMAKE_SOURCE_DIR}/include)
endforeach (DEMO)

