find_path(Julia_INCLUDE_DIRS
	NAMES julia.h
	PATH_SUFFIXES julia
	HINTS /usr/local/include /usr/include /opt/local/include
	DOC "Julia include directories"
)

list(APPEND Julia_INCLUDE_DIRS "${Julia_INCLUDE_DIRS}/..")


find_library(Julia_LIBRARIES
	NAMES julia
	HINTS /usr/local/lib64 /usr/local/lib /usr/lib64 /usr/lib /opt/local/lib /usr/lib32 /usr/local/lib32
	DOC "Julia libraries"
)

if(Julia_INCLUDE_DIRS AND Julia_LIBRARIES)
	set(Julia_FOUND TRUE)

	message("Julia include directories: ${Julia_INCLUDE_DIRS}")
	message("Julia library: ${Julia_LIBRARIES}")
else()
	set(Julia_FOUND FALSE)

	message("Error: Julia could not be found!")
endif()
