#
# finds the minuit libs
# @author Tobias Weber <tobias.weber@tum.de>
# @license GPLv2 or GPLv3
#

find_path(Minuit2_INCLUDE_DIRS
	NAMES MnMigrad.h
	PATH_SUFFIXES root Minuit2 Minuit root/Minuit2 root/Minuit Minuit2/Minuit2
	HINTS /usr/local/include/ /usr/include/ /opt/local/include
	DOC "Root/Minuit2 include directories"
)

# also include root base dir
list(APPEND Minuit2_INCLUDE_DIRS "${Minuit2_INCLUDE_DIRS}/..")


find_library(Minuit2_LIBRARIES
	NAMES Minuit2
	HINTS /usr/local/lib64 /usr/local/lib64/root /usr/local/lib /usr/local/lib/root /usr/lib64 /usr/lib64/root /usr/lib /usr/lib/root /opt/local/lib /opt/local/lib/root /usr/lib32 /usr/lib32/root
	DOC "Minuit2 library"
)

if(Minuit2_INCLUDE_DIRS AND Minuit2_LIBRARIES)
	set(Minuit2_FOUND TRUE)

	message("Minuit include directories: ${Minuit2_INCLUDE_DIRS}")
	message("Minuit library: ${Minuit2_LIBRARIES}")
else()
	set(Minuit2_FOUND FALSE)

	message("Error: Minuit2 could not be found!")
endif()
