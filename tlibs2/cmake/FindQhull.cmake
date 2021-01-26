#
# finds the qhull libs
# @author Tobias Weber <tweber@ill.fr>
# @date 13-aug-2020
# @license GPLv3
#

find_path(Qhull_INCLUDE_DIRS
	NAMES Qhull.h
	PATH_SUFFIXES libqhullcpp libqhull_r qhull
	HINTS /usr/local/include/ /usr/include/ /opt/local/include
	DOC "Qhull include directories"
)


find_library(Qhull_c_LIBRARY
	NAMES qhull_r
	HINTS /usr/local/lib64 /usr/local/lib /usr/lib64 /usr/lib /opt/local/lib
	DOC "Qhull C library"
)


find_library(Qhull_cpp_LIBRARY
	NAMES qhullcpp
	HINTS /usr/local/lib64 /usr/local/lib /usr/lib64 /usr/lib /opt/local/lib
	DOC "Qhull C++ library"
)


if(Qhull_INCLUDE_DIRS AND Qhull_c_LIBRARY AND Qhull_cpp_LIBRARY)
	set(Qhull_FOUND TRUE)

	list(APPEND Qhull_LIBRARIES "${Qhull_c_LIBRARY}")
	list(APPEND Qhull_LIBRARIES "${Qhull_cpp_LIBRARY}")

	message("Qhull include directories: ${Qhull_INCLUDE_DIRS}")
	message("Qhull libraries: ${Qhull_LIBRARIES}")
else()
	set(Qhull_FOUND FALSE)

	message("Error: Qhull could not be found!")
endif()
