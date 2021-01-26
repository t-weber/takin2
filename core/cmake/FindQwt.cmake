find_path(QWT_INCLUDE_DIRS
	NAMES qwt.h
	PATH_SUFFIXES qwt6-qt5 qwt-qt5 qt5/qwt qt5/qwt/qwt qt5/qwt6 qt5 qwt6 qwt
	HINTS /usr/local/include /usr/include /opt/local/include /usr/local/Cellar/qwt/*/lib/qwt.framework/Versions/6/Headers/
	DOC "Qwt include directories"
)

list(APPEND QWT_INCLUDE_DIRS "${QWT_INCLUDE_DIRS}/..")


find_library(QWT_LIBRARIES
	NAMES qwt6-qt5 qwt-qt5 qwt6 qwt
	HINTS /usr/local/lib64 /usr/local/lib /usr/lib64 /usr/lib /opt/local/lib /usr/lib32 /usr/local/lib32
	DOC "Qwt libraries"
)

if(QWT_INCLUDE_DIRS AND QWT_LIBRARIES)
	set(QWT_FOUND TRUE)

	message("Qwt include directories: ${QWT_INCLUDE_DIRS}")
	message("Qwt library: ${QWT_LIBRARIES}")
else()
	set(QWT_FOUND FALSE)

	message("Error: Qwt could not be found!")
endif()
