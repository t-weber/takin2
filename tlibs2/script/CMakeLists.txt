#
# @author Tobias Weber <tweber@ill.fr>
# @date 11-apr-2020
# @license: see 'LICENSE' file
# @desc Forked on 18/July/2020 from my privately developed "matrix_calc" project (https://github.com/t-weber/matrix_calc).
#

project(mcalc)
cmake_minimum_required(VERSION 3.0)

list(APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}" "${PROJECT_SOURCE_DIR}/../cmake")

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_VERBOSE_MAKEFILE ON)

find_package(FLEX 2 REQUIRED)
find_package(BISON 3 REQUIRED)
find_package(Boost REQUIRED COMPONENTS program_options)
find_package(Lapacke REQUIRED)

add_compile_options(-Wall -Wextra)
add_compile_options(-std=c++20)

add_definitions(${Boost_CXX_FLAGS})
add_definitions(-DUSE_LAPACK)


include_directories(
	"${PROJECT_SOURCE_DIR}"
	"${PROJECT_SOURCE_DIR}/src" "${PROJECT_SOURCE_DIR}/libs"
	"${PROJECT_SOURCE_DIR}/.."
	"${FLEX_INCLUDE_DIRS}" "${Boost_INCLUDE_DIRS}" "${Boost_INCLUDE_DIRS}/.."
	"${CMAKE_CURRENT_BINARY_DIR}"
	"${Lapacke_INCLUDE_DIRS}"
)


BISON_TARGET(parser_impl
	src/parser.y ${CMAKE_CURRENT_BINARY_DIR}/parser_impl.cpp
	DEFINES_FILE ${CMAKE_CURRENT_BINARY_DIR}/parser_defs.h
	VERBOSE REPORT_FILE ${CMAKE_CURRENT_BINARY_DIR}/parser_report.txt
)


FLEX_TARGET(lexer_impl
	src/lexer.l ${CMAKE_CURRENT_BINARY_DIR}/lexer_impl.cpp
	DEFINES_FILE ${CMAKE_CURRENT_BINARY_DIR}/lexer_impl.h
)

ADD_FLEX_BISON_DEPENDENCY(lexer_impl parser_impl)



message("Flex outputs: ${FLEX_lexer_impl_OUTPUTS}")
message("Bison output source: ${BISON_parser_impl_OUTPUT_SOURCE}")
message("Bison output header: ${BISON_parser_impl_OUTPUT_HEADER}")



add_executable(mcalc
	src/parser.cpp src/parser.h src/ast.h src/sym.h
	src/llasm.cpp src/llasm_ops.cpp src/llasm_var.cpp
	src/llasm_arr.cpp src/llasm_func.cpp src/llasm.h
	src/printast.cpp src/printast.h
	${FLEX_lexer_impl_OUTPUTS}
	${BISON_parser_impl_OUTPUT_SOURCE} ${BISON_parser_impl_OUTPUT_HEADER}
	../libs/log.cpp ../libs/log.h
)

#add_dependencies(parser parser_impl lexer_impl)
target_link_libraries(mcalc ${Boost_LIBRARIES})




add_library(mcalc_rt SHARED
	src/runtime.cpp
)

target_link_libraries(mcalc_rt ${Lapacke_LIBRARIES})