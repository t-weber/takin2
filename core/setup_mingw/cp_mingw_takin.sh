#!/bin/bash
#
# creates a distro for mingw
# @author Tobias Weber <tobias.weber@tum.de>
# @date 2016-2020
# @license GPLv2
#


# installation directory
INSTDIR="$1"

if [ "${INSTDIR}" = "" ]; then
	INSTDIR=~/.wine/drive_c/takin
fi

mkdir -p ${INSTDIR}


# main programs
cp -v bin/*.exe		${INSTDIR}/


# info files
cp -v *.txt			${INSTDIR}/
cp -rv 3rdparty_licenses/	${INSTDIR}/


# examples
cp -rv examples 		${INSTDIR}/
cp -rv data/samples 		${INSTDIR}/
cp -rv data/instruments 	${INSTDIR}/
cp -rv data/demos 		${INSTDIR}/


# resources
mkdir ${INSTDIR}/res
cp -rv res/* ${INSTDIR}/res/
cp -rv doc/*.html ${INSTDIR}/res/doc/
gunzip -v ${INSTDIR}/res/data/*



# libraries
cp -v /usr/x86_64-w64-mingw32/sys-root/mingw/bin/libstdc++-6.dll	${INSTDIR}/
cp -v /usr/x86_64-w64-mingw32/sys-root/mingw/bin/libwinpthread-1.dll	${INSTDIR}/
cp -v /usr/x86_64-w64-mingw32/sys-root/mingw/bin/libpng16-16.dll	${INSTDIR}/
cp -v /usr/x86_64-w64-mingw32/sys-root/mingw/bin/libgcc_s_sjlj-1.dll	${INSTDIR}/
cp -v /usr/x86_64-w64-mingw32/sys-root/mingw/bin/libgcc_s_seh-1.dll	${INSTDIR}/
cp -v /usr/x86_64-w64-mingw32/sys-root/mingw/bin/libbz2-1.dll		${INSTDIR}/
cp -v /usr/x86_64-w64-mingw32/sys-root/mingw/bin/zlib1.dll		${INSTDIR}/
cp -v /usr/x86_64-w64-mingw32/sys-root/mingw/bin/iconv.dll		${INSTDIR}/
cp -v /usr/x86_64-w64-mingw32/sys-root/mingw/bin/libpcre2-16-0.dll	${INSTDIR}/
cp -v /usr/x86_64-w64-mingw32/sys-root/mingw/bin/libharfbuzz-0.dll	${INSTDIR}/
cp -v /usr/x86_64-w64-mingw32/sys-root/mingw/bin/libglib-2.0-0.dll	${INSTDIR}/
cp -v /usr/x86_64-w64-mingw32/sys-root/mingw/bin/libintl-8.dll		${INSTDIR}/
cp -v /usr/x86_64-w64-mingw32/sys-root/mingw/bin/libpcre-1.dll		${INSTDIR}/
cp -v /usr/x86_64-w64-mingw32/sys-root/mingw/bin/libssp-0.dll		${INSTDIR}/
#cp -v /usr/x86_64-w64-mingw32/sys-root/mingw/lib/libMinuit2.dll	${INSTDIR}/

cp -v /usr/x86_64-w64-mingw32/sys-root/mingw/bin/libboost_regex-x64.dll			${INSTDIR}/
cp -v /usr/x86_64-w64-mingw32/sys-root/mingw/bin/libboost_system-x64.dll		${INSTDIR}/
cp -v /usr/x86_64-w64-mingw32/sys-root/mingw/bin/libboost_iostreams-x64.dll		${INSTDIR}/
cp -v /usr/x86_64-w64-mingw32/sys-root/mingw/bin/libboost_filesystem-x64.dll		${INSTDIR}/
cp -v /usr/x86_64-w64-mingw32/sys-root/mingw/bin/libboost_program_options-x64.dll	${INSTDIR}/

#cp -v /usr/x86_64-w64-mingw32/sys-root/mingw/bin/libboost_python39.dll	${INSTDIR}/
#cp -v /usr/x86_64-w64-mingw32/sys-root/mingw/bin/libpython3.9.dll	${INSTDIR}/
#cp -v /usr/x86_64-w64-mingw32/sys-root/mingw/bin/libcrypto-1_1-x64.dll	${INSTDIR}/
#cp -v /usr/x86_64-w64-mingw32/sys-root/mingw/bin/libffi-6.dll		${INSTDIR}/

cp -v /usr/x86_64-w64-mingw32/sys-root/mingw/bin/Qt5Core.dll		${INSTDIR}/
cp -v /usr/x86_64-w64-mingw32/sys-root/mingw/bin/Qt5Gui.dll		${INSTDIR}/
cp -v /usr/x86_64-w64-mingw32/sys-root/mingw/bin/Qt5Widgets.dll		${INSTDIR}/
cp -v /usr/x86_64-w64-mingw32/sys-root/mingw/bin/Qt5OpenGL.dll		${INSTDIR}/
cp -v /usr/x86_64-w64-mingw32/sys-root/mingw/bin/Qt5Svg.dll		${INSTDIR}/
cp -v /usr/x86_64-w64-mingw32/sys-root/mingw/bin/Qt5Xml.dll		${INSTDIR}/
cp -v /usr/x86_64-w64-mingw32/sys-root/mingw/bin/Qt5PrintSupport.dll	${INSTDIR}/
cp -v /usr/x86_64-w64-mingw32/sys-root/mingw/bin/qwt-qt5.dll		${INSTDIR}/
cp -v /usr/x86_64-w64-mingw32/sys-root/mingw/bin/libfreetype-6.dll	${INSTDIR}/

#cp -v /usr/x86_64-w64-mingw32/sys-root/mingw/bin/liblapack.dll		${INSTDIR}/
#cp -v /usr/x86_64-w64-mingw32/sys-root/mingw/bin/liblapacke.dll	${INSTDIR}/
#cp -v /usr/x86_64-w64-mingw32/sys-root/mingw/bin/libblas.dll		${INSTDIR}/
#cp -v /usr/x86_64-w64-mingw32/sys-root/mingw/bin/libquadmath-0.dll	${INSTDIR}/
#cp -v /usr/x86_64-w64-mingw32/sys-root/mingw/bin/libgfortran-5.dll	${INSTDIR}/


# qt plugins
mkdir -p ${INSTDIR}/lib/plugins/platforms/
mkdir -p ${INSTDIR}/lib/plugins/iconengines/
mkdir -p ${INSTDIR}/lib/plugins/imageformats

cp -v /usr/x86_64-w64-mingw32/sys-root/mingw/lib/qt5/plugins/platforms/*.dll		${INSTDIR}/lib/plugins/platforms/
cp -v /usr/x86_64-w64-mingw32/sys-root/mingw/lib/qt5/plugins/iconengines/qsvgicon.dll	${INSTDIR}/lib/plugins/iconengines/
cp -v /usr/x86_64-w64-mingw32/sys-root/mingw/lib/qt5/plugins/imageformats/qsvg.dll	${INSTDIR}/lib/plugins/imageformats/



# stripping
strip -v ${INSTDIR}/*.exe
strip -v ${INSTDIR}/*.dll
strip -v ${INSTDIR}/lib/plugins/platforms/*.dll
strip -v ${INSTDIR}/lib/plugins/iconengines/*.dll
strip -v ${INSTDIR}/lib/plugins/imageformats/*.dll


echo -e "[Paths]\nPlugins = lib/plugins\n" > ${INSTDIR}/qt.conf
