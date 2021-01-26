#
# finds the lapacke libs
# @author Tobias Weber <tweber@ill.fr>
# @date 8-aug-2020
# @license GPLv3
#

find_path(Lapacke_INCLUDE_DIRS
	NAMES lapacke.h
	PATH_SUFFIXES lapacke lapacke/include
	HINTS /usr/local/include/ /usr/include/ /opt/local/include /usr/local/opt/lapack/include /usr/local/Cellar/lapack/*/include
	DOC "Lapacke include directories"
)


find_library(Lapacke_LIBRARIES
	NAMES lapacke
	HINTS /usr/local/lib64 /usr/local/lib /usr/lib64 /usr/lib /opt/local/lib /usr/local/opt/lapack/lib /usr/local/Cellar/lapack/*/lib
	DOC "Lapacke library"
)


if(Lapacke_INCLUDE_DIRS AND Lapacke_LIBRARIES)
	set(Lapacke_FOUND TRUE)

	message("Lapacke include directories: ${Lapacke_INCLUDE_DIRS}")
	message("Lapacke library: ${Lapacke_LIBRARIES}")
else()
	set(Lapacke_FOUND FALSE)

	message("Error: Lapacke could not be found!")
endif()
