#
# tlibs2 -- tests
# @author Tobias Weber <tweber@ill.fr>
# @date 8-jun-2020
# @license see 'LICENSE' file
#

project(tlibs2_tests)

cmake_minimum_required(VERSION 3.0)
cmake_policy(SET CMP0078 NEW)
cmake_policy(SET CMP0086 NEW)
enable_testing()

list(APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}" "${PROJECT_SOURCE_DIR}/../cmake")


message("Build type: ${CMAKE_BUILD_TYPE}")

if(NOT "${CMAKE_BUILD_TYPE}" STREQUAL "Release")
	set(CMAKE_VERBOSE_MAKEFILE TRUE)
endif()


set(CMAKE_CXX_STANDARD 20)
add_compile_options(-std=c++20)
add_compile_options(-Wall -Wextra)


# -----------------------------------------------------------------------------
# Boost
add_definitions(-DBOOST_SYSTEM_NO_DEPRECATED)
find_package(Boost REQUIRED)
add_definitions(${Boost_CXX_FLAGS})

message("Using Boost version ${Boost_MAJOR_VERSION}.${Boost_MINOR_VERSION}.${Boost_SUBMINOR_VERSION}.")
# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
# Lapacke
find_package(Lapacke REQUIRED)
add_definitions(-DUSE_LAPACK)
# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
# further packages
find_package(Minuit2 REQUIRED)
find_package(Threads REQUIRED)
find_package(OpenMP REQUIRED)
# -----------------------------------------------------------------------------


include_directories(
	"${PROJECT_SOURCE_DIR}" "${PROJECT_SOURCE_DIR}/.."
	"${Minuit2_INCLUDE_DIRS}"
	"${Lapacke_INCLUDE_DIRS}"
)


# -----------------------------------------------------------------------------
# build individual test cases
add_executable(mat0 mat0.cpp)
add_executable(mat1 mat1.cpp)
add_executable(stat1 stat1.cpp)
add_executable(quat1 quat1.cpp)
add_executable(expr expr.cpp ../libs/log.cpp)
add_executable(cryst cryst.cpp ../libs/log.cpp)
add_executable(leastsq leastsq.cpp ../libs/log.cpp)
add_executable(eig eig.cpp)
add_executable(ode ode.cpp)
add_executable(min0 min0.cpp ../libs/log.cpp)
add_executable(fit0 fit0.cpp ../libs/log.cpp)
add_executable(fit1 fit1.cpp ../libs/log.cpp)
add_executable(cov cov.cpp)

target_link_libraries(expr Threads::Threads)
target_link_libraries(mat0 ${Lapacke_LIBRARIES})
target_link_libraries(cryst ${Lapacke_LIBRARIES})
target_link_libraries(leastsq ${Lapacke_LIBRARIES})
target_link_libraries(eig ${Lapacke_LIBRARIES})
target_link_libraries(ode ${Lapacke_LIBRARIES})
target_link_libraries(min0 ${Minuit2_LIBRARIES} Threads::Threads OpenMP::OpenMP_CXX)
target_link_libraries(fit0 ${Minuit2_LIBRARIES} Threads::Threads OpenMP::OpenMP_CXX)
target_link_libraries(fit1 ${Minuit2_LIBRARIES} OpenMP::OpenMP_CXX)
# -----------------------------------------------------------------------------


# -----------------------------------------------------------------------------
# run individual test cases
add_test(expr expr)
add_test(cryst cryst)
add_test(mat0 mat0)
add_test(mat1 mat1)
add_test(stat1 stat1)
add_test(quat1 quat1)
add_test(lestsq leastsq)
add_test(eig eig)
add_test(ode ode)
add_test(min0 min0)
add_test(fit0 fit0)
add_test(fit1 fit1)
add_test(cov cov)
# -----------------------------------------------------------------------------
