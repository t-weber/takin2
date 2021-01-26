#
# random distributions
# @author Tobias Weber <tobias.weber@tum.de>
# @date sep-2016
# @license GPLv2 or GPLv3
#

project(distr)
cmake_minimum_required(VERSION 3.0)


set(CMAKE_VERBOSE_MAKEFILE TRUE)
#set(CMAKE_BUILD_TYPE Release)
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_AUTOMOC ON)
set(CMAKE_AUTOUIC ON)

list(APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}" ".")


find_package(Boost REQUIRED)
find_package(Threads REQUIRED)
find_package(Qt5Core REQUIRED)
find_package(Qt5Gui REQUIRED)
find_package(Qt5Widgets REQUIRED)

add_definitions(-I${CMAKE_CURRENT_BINARY_DIR})
add_definitions(-DINSTALL_PREFIX="${CMAKE_INSTALL_PREFIX}")
add_definitions(${Boost_CXX_FLAGS})
add_definitions(-std=c++14)


include_directories("${PROJECT_SOURCE_DIR}" "." 
	"${PROJECT_SOURCE_DIR}/.." "${PROJECT_SOURCE_DIR}/../..")
include_directories("${Boost_INCLUDE_DIRS}/..")
include_directories("/usr/include/qwt")


add_executable(distr 
	../../log/log.cpp
	distr_main.cpp distr.cpp distr.h)

target_link_libraries(distr 
	${Boost_LIBRARIES} Threads::Threads
	Qt5::Core Qt5::Widgets Qt5::Gui
	qwt-qt5)
