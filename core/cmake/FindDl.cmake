find_library(Dl_LIBRARIES
	NAMES libdl.so libdl.so.2
	HINTS /usr/local/lib64 /usr/local/lib /usr/lib64 /usr/lib /opt/local/lib /lib/x86_64-linux-gnu #/usr/lib32 /usr/local/lib32
	DOC "DL library"
)

if(Dl_LIBRARIES)
	set(Dl_FOUND TRUE)
	message("DL library: ${Dl_LIBRARIES}")
else()
	set(Dl_FOUND FALSE)
	message("Error: libdl could not be found!")
endif()
