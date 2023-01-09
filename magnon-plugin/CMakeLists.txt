#
# module for magnon dynamics
# @author Tobias Weber <tweber@ill.fr>
# @date 28-nov-2022
# @license see 'LICENSE' file
#
# ----------------------------------------------------------------------------
# magnonmod (part of the Takin software suite)
# Copyright (C) 2022  Tobias WEBER (Institut Laue-Langevin (ILL),
#                     Grenoble, France).
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, version 2 of the License.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
# ----------------------------------------------------------------------------
#

cmake_minimum_required(VERSION 3.0)
project(magnonmod)
enable_language(CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_POSITION_INDEPENDENT_CODE TRUE)
add_compile_options(-std=c++20)
add_compile_options(-Wall -Wextra -Weffc++)
add_compile_options(-fvisibility=hidden)

list(APPEND CMAKE_MODULE_PATH
	"${PROJECT_SOURCE_DIR}"
	"${PROJECT_SOURCE_DIR}/ext/tlibs2/cmake"
)


if(CMAKE_BUILD_TYPE STREQUAL "Debug")
	add_compile_options(-g -ggdb)
	set(CMAKE_VERBOSE_MAKEFILE TRUE)

	message("Debug building for ${CMAKE_SYSTEM_NAME} systems.")

elseif(CMAKE_BUILD_TYPE STREQUAL "Release")
	add_compile_options("-DNDEBUG")
	add_compile_options("-Wno-#pragma-messages")

	if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
		# see: https://gcc.gnu.org/onlinedocs/gcc/Developer-Options.html
		add_compile_options(-fcompare-debug-second)
	endif()

	message("Release building for ${CMAKE_SYSTEM_NAME} systems.")
endif()


if(CMAKE_SYSTEM_NAME STREQUAL "Windows")
	# pass linker --subsystem option
	add_compile_options(-Wl,--subsystem,windows)
elseif(CMAKE_SYSTEM_NAME STREQUAL "Darwin")
	add_compile_options(-mmacosx-version-min=10.15)
endif()


find_package(Boost REQUIRED)
add_definitions(${Boost_CXX_FLAGS})

find_package(Lapacke REQUIRED)
add_definitions(-DUSE_LAPACK=1)
include_directories("${Lapacke_INCLUDE_DIRS}")

include_directories("${PROJECT_SOURCE_DIR}")

include_directories(
	SYSTEM "${Boost_INCLUDE_DIRS}" "${Boost_INCLUDE_DIRS}/.."
	"../.." "ext" "ext/takin"
)


add_library(magnonmod SHARED
	src/magnonmod.cpp src/magnonmod.h
	ext/tlibs2/libs/magdyn.h
	ext/tlibs2/libs/maths.h

	ext/tlibs2/libs/log.cpp ext/tlibs/log/log.cpp
	ext/takin/tools/monteconvo/sqwbase.cpp
)

if(CMAKE_SYSTEM_NAME STREQUAL "Darwin")
	target_link_options(magnonmod
		# create an __info_plist section in the binary
		PRIVATE LINKER:-sectcreate,__TEXT,__info_plist,${PROJECT_SOURCE_DIR}/ext/takin2/setup_mac/InfoMagDyn.plist
	)
endif()

target_link_libraries(magnonmod
	${Boost_LIBRARIES}
	${Lapacke_LIBRARIES}
)
