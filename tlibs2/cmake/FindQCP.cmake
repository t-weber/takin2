#
# finds the qcp libs
# @author Tobias Weber <tweber@ill.fr>
# @date 16-sep-2020
# @license GPLv3
#

find_path(QCP_INCLUDE_DIRS
	NAMES qcustomplot.h
	HINTS /usr/local/include/ /usr/include/ /opt/local/include
	DOC "QCP include directories"
)


find_library(QCP_LIBRARIES
	NAMES qcustomplot
	HINTS /usr/local/lib64 /usr/local/lib /usr/lib/x86_64-linux-gnu /usr/lib64 /usr/lib /opt/local/lib
	DOC "QCP library"
)


if(QCP_INCLUDE_DIRS AND QCP_LIBRARIES)
	set(QCP_FOUND TRUE)

	message("QCP include directories: ${QCP_INCLUDE_DIRS}")
	message("QCP libraries: ${QCP_LIBRARIES}")
else()
	set(QCP_FOUND FALSE)

	message("Error: QCP could not be found!")
endif()
