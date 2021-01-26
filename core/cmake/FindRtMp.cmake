find_library(Rt_LIBRARIES
	NAMES librt.so
	HINTS /usr/local/lib64 /usr/local/lib /usr/lib64 /usr/lib /opt/local/lib /lib/x86_64-linux-gnu #/usr/lib32 /usr/local/lib32
	DOC "RT library"
)

if(Rt_LIBRARIES)
	set(Rt_FOUND TRUE)
	message("RT library: ${Rt_LIBRARIES}")
else()
	set(Rt_FOUND FALSE)
	message("Error: librt could not be found!")
endif()



find_library(Mp_LIBRARIES
	NAMES libgomp.so libgomp.so.1
	HINTS /usr/local/lib64 /usr/local/lib /usr/lib64 /usr/lib /opt/local/lib /lib/x86_64-linux-gnu #/usr/lib32 /usr/local/lib32
	DOC "MP library"
)

if(Mp_LIBRARIES)
	set(Mp_FOUND TRUE)
	message("MP library: ${Mp_LIBRARIES}")
else()
	set(Mp_FOUND FALSE)
	message("Error: libgomp could not be found!")
endif()
