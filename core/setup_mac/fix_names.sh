#!/bin/bash
#
# @date 22-mar-17
# @author Tobias Weber <tobias.weber@tum.de>
# @license GPLv2
#
# makes the framework libraries locally usable
#

PRG="takin.app"

TOOL=install_name_tool
STRIP=strip

QT_VER=$(ls /usr/local/Cellar/qt)
PY_VER=3.9

echo -e "Qt version: ${QT_VER}"
echo -e "Py version: ${PY_VER}"


# files whose linkage is to be changed
declare -a filestochange=(
	"${PRG}/Contents/Frameworks/QtCore.framework/Versions/5/QtCore"
	"${PRG}/Contents/Frameworks/QtWidgets.framework/Versions/5/QtWidgets"
	"${PRG}/Contents/Frameworks/QtGui.framework/Versions/5/QtGui"
	"${PRG}/Contents/Frameworks/QtConcurrent.framework/Versions/5/QtConcurrent"
	"${PRG}/Contents/Frameworks/QtOpenGL.framework/Versions/5/QtOpenGL"
	"${PRG}/Contents/Frameworks/QtPrintSupport.framework/Versions/5/QtPrintSupport"
	"${PRG}/Contents/Frameworks/QtDBus.framework/Versions/5/QtDBus"
	"${PRG}/Contents/Frameworks/QtSvg.framework/Versions/5/QtSvg"
	"${PRG}/Contents/Frameworks/QtXml.framework/Versions/5/QtXml"
	"${PRG}/Contents/Frameworks/QtXmlPatterns.framework/Versions/5/QtXmlPatterns"
	"${PRG}/Contents/Frameworks/qwt.framework/Versions/6/qwt"
	"${PRG}/Contents/MacOS/takin"
	"${PRG}/Contents/MacOS/takin_cif2xml"
	"${PRG}/Contents/MacOS/takin_findsg"
	"${PRG}/Contents/MacOS/takin_pol"
	"${PRG}/Contents/MacOS/takin_structfact"
	"${PRG}/Contents/MacOS/takin_magstructfact"
	"${PRG}/Contents/MacOS/takin_scanbrowser"
	"${PRG}/Contents/MacOS/takin_magsgbrowser"
	"${PRG}/Contents/MacOS/takin_moldyn"
	"${PRG}/Contents/MacOS/takinmod_py"
	"${PRG}/Contents/MacOS/takinmod_jl"
	"${PRG}/Contents/MacOS/takin_convofit"
	"${PRG}/Contents/MacOS/takin_convoseries"
	"${PRG}/Contents/MacOS/takin_polextract"
)



# original symbols
declare -a changefrom=(
	"/usr/local/Cellar/qt/${QT_VER}/lib/QtCore.framework/Versions/5/QtCore"
	"/usr/local/Cellar/qt5/${QT_VER}/lib/QtCore.framework/Versions/5/QtCore"
	"/usr/local/opt/qt/lib/QtCore.framework/Versions/5/QtCore"
	"/usr/local/opt/qt5/lib/QtCore.framework/Versions/5/QtCore"
	"/usr/local/Cellar/qt/${QT_VER}/lib/QtGui.framework/Versions/5/QtGui"
	"/usr/local/Cellar/qt5/${QT_VER}/lib/QtGui.framework/Versions/5/QtGui"
	"/usr/local/opt/qt/lib/QtGui.framework/Versions/5/QtGui"
	"/usr/local/opt/qt5/lib/QtGui.framework/Versions/5/QtGui"
	"/usr/local/Cellar/qt/${QT_VER}/lib/QtWidgets.framework/Versions/5/QtWidgets"
	"/usr/local/Cellar/qt5/${QT_VER}/lib/QtWidgets.framework/Versions/5/QtWidgets"
	"/usr/local/opt/qt/lib/QtWidgets.framework/Versions/5/QtWidgets"
	"/usr/local/opt/qt5/lib/QtWidgets.framework/Versions/5/QtWidgets"
	"/usr/local/Cellar/qt/${QT_VER}/lib/QtOpenGL.framework/Versions/5/QtOpenGL"
	"/usr/local/Cellar/qt5/${QT_VER}/lib/QtOpenGL.framework/Versions/5/QtOpenGL"
	"/usr/local/opt/qt/lib/QtOpenGL.framework/Versions/5/QtOpenGL"
	"/usr/local/opt/qt5/lib/QtOpenGL.framework/Versions/5/QtOpenGL"
	"/usr/local/Cellar/qt/${QT_VER}/lib/QtConcurrent.framework/Versions/5/QtConcurrent"
	"/usr/local/Cellar/qt5/${QT_VER}/lib/QtConcurrent.framework/Versions/5/QtConcurrent"
	"/usr/local/opt/qt/lib/QtConcurrent.framework/Versions/5/QtConcurrent"
	"/usr/local/opt/qt5/lib/QtConcurrent.framework/Versions/5/QtConcurrent"
	"/usr/local/Cellar/qt/${QT_VER}/lib/QtXml.framework/Versions/5/QtXml"
	"/usr/local/Cellar/qt5/${QT_VER}/lib/QtXml.framework/Versions/5/QtXml"
	"/usr/local/opt/qt/lib/QtXml.framework/Versions/5/QtXml"
	"/usr/local/opt/qt5/lib/QtXml.framework/Versions/5/QtXml"
	"/usr/local/Cellar/qt/${QT_VER}/lib/QtXmlPatterns.framework/Versions/5/QtXmlPatterns"
	"/usr/local/Cellar/qt5/${QT_VER}/lib/QtXmlPatterns.framework/Versions/5/QtXmlPatterns"
	"/usr/local/opt/qt/lib/QtXmlPatterns.framework/Versions/5/QtXmlPatterns"
	"/usr/local/opt/qt5/lib/QtXmlPatterns.framework/Versions/5/QtXmlPatterns"
	"/usr/local/Cellar/qt/${QT_VER}/lib/QtSvg.framework/Versions/5/QtSvg"
	"/usr/local/Cellar/qt5/${QT_VER}/lib/QtSvg.framework/Versions/5/QtSvg"
	"/usr/local/opt/qt/lib/QtSvg.framework/Versions/5/QtSvg"
	"/usr/local/opt/qt5/lib/QtSvg.framework/Versions/5/QtSvg"
	"/usr/local/Cellar/qt/${QT_VER}/lib/QtPrintSupport.framework/Versions/5/QtPrintSupport"
	"/usr/local/Cellar/qt5/${QT_VER}/lib/QtPrintSupport.framework/Versions/5/QtPrintSupport"
	"/usr/local/opt/qt/lib/QtPrintSupport.framework/Versions/5/QtPrintSupport"
	"/usr/local/opt/qt5/lib/QtPrintSupport.framework/Versions/5/QtPrintSupport"
	"/usr/local/Cellar/qt/${QT_VER}/lib/QtDBus.framework/Versions/5/QtDBus"
	"/usr/local/Cellar/qt5/${QT_VER}/lib/QtDBus.framework/Versions/5/QtDBus"
	"/usr/local/opt/qt/lib/QtDBus.framework/Versions/5/QtDBus"
	"/usr/local/opt/qt5/lib/QtDBus.framework/Versions/5/QtDBus"
	"/usr/local/opt/qwt/lib/qwt.framework/Versions/6/qwt"
	"/usr/local/opt/python@${PY_VER}/Frameworks/Python.framework/Versions/${PY_VER}/Python"
	"/usr/local/opt/minuit2/lib/libMinuit2.0.dylib"
	"/usr/local/opt/boost/lib/libboost_system-mt.dylib"
	"/usr/local/opt/boost/lib/libboost_filesystem-mt.dylib"
	"/usr/local/opt/boost/lib/libboost_iostreams-mt.dylib"
	"/usr/local/opt/boost/lib/libboost_program_options-mt.dylib"
	"/usr/local/opt/boost-python3/lib/libboost_python39-mt.dylib"
	"/usr/local/opt/freetype/lib/libfreetype.6.dylib"
	"/usr/local/opt/libpng/lib/libpng16.16.dylib"
	"/usr/local/opt/libjpeg/lib/libjpeg.9.dylib"
	"/usr/local/opt/jpeg/lib/libjpeg.9.dylib"
	"/usr/local/opt/libtiff/lib/libtiff.5.dylib"
	"/usr/local/opt/openblas/lib/libopenblas.0.dylib"
	"/usr/local/opt/gcc/lib/gcc/10/libgfortran.5.dylib"
	"/usr/local/opt/gcc/lib/gcc/10/libgomp.1.dylib"
	"/usr/local/opt/gcc/lib/gcc/10/libquadmath.0.dylib"
	"/usr/local/Cellar/gcc/10.2.0/lib/gcc/10/libquadmath.0.dylib"
	"/usr/local/opt/gcc/lib/gcc/10/libgcc_s.1.dylib"
	"/usr/local/lib/gcc/10/libgcc_s.1.dylib"
)


# symbols to change into
declare -a changeto=(
	"@executable_path/../Frameworks/QtCore.framework/Versions/5/QtCore"
	"@executable_path/../Frameworks/QtCore.framework/Versions/5/QtCore"
	"@executable_path/../Frameworks/QtCore.framework/Versions/5/QtCore"
	"@executable_path/../Frameworks/QtCore.framework/Versions/5/QtCore"
	"@executable_path/../Frameworks/QtGui.framework/Versions/5/QtGui"
	"@executable_path/../Frameworks/QtGui.framework/Versions/5/QtGui"
	"@executable_path/../Frameworks/QtGui.framework/Versions/5/QtGui"
	"@executable_path/../Frameworks/QtGui.framework/Versions/5/QtGui"
	"@executable_path/../Frameworks/QtWidgets.framework/Versions/5/QtWidgets"
	"@executable_path/../Frameworks/QtWidgets.framework/Versions/5/QtWidgets"
	"@executable_path/../Frameworks/QtWidgets.framework/Versions/5/QtWidgets"
	"@executable_path/../Frameworks/QtWidgets.framework/Versions/5/QtWidgets"
	"@executable_path/../Frameworks/QtOpenGL.framework/Versions/5/QtOpenGL"
	"@executable_path/../Frameworks/QtOpenGL.framework/Versions/5/QtOpenGL"
	"@executable_path/../Frameworks/QtOpenGL.framework/Versions/5/QtOpenGL"
	"@executable_path/../Frameworks/QtOpenGL.framework/Versions/5/QtOpenGL"
	"@executable_path/../Frameworks/QtConcurrent.framework/Versions/5/QtConcurrent"
	"@executable_path/../Frameworks/QtConcurrent.framework/Versions/5/QtConcurrent"
	"@executable_path/../Frameworks/QtConcurrent.framework/Versions/5/QtConcurrent"
	"@executable_path/../Frameworks/QtConcurrent.framework/Versions/5/QtConcurrent"
	"@executable_path/../Frameworks/QtXml.framework/Versions/5/QtXml"
	"@executable_path/../Frameworks/QtXml.framework/Versions/5/QtXml"
	"@executable_path/../Frameworks/QtXml.framework/Versions/5/QtXml"
	"@executable_path/../Frameworks/QtXml.framework/Versions/5/QtXml"
	"@executable_path/../Frameworks/QtXmlPatterns.framework/Versions/5/QtXmlPatterns"
	"@executable_path/../Frameworks/QtXmlPatterns.framework/Versions/5/QtXmlPatterns"
	"@executable_path/../Frameworks/QtXmlPatterns.framework/Versions/5/QtXmlPatterns"
	"@executable_path/../Frameworks/QtXmlPatterns.framework/Versions/5/QtXmlPatterns"
	"@executable_path/../Frameworks/QtSvg.framework/Versions/5/QtSvg"
	"@executable_path/../Frameworks/QtSvg.framework/Versions/5/QtSvg"
	"@executable_path/../Frameworks/QtSvg.framework/Versions/5/QtSvg"
	"@executable_path/../Frameworks/QtSvg.framework/Versions/5/QtSvg"
	"@executable_path/../Frameworks/QtPrintSupport.framework/Versions/5/QtPrintSupport"
	"@executable_path/../Frameworks/QtPrintSupport.framework/Versions/5/QtPrintSupport"
	"@executable_path/../Frameworks/QtPrintSupport.framework/Versions/5/QtPrintSupport"
	"@executable_path/../Frameworks/QtPrintSupport.framework/Versions/5/QtPrintSupport"
	"@executable_path/../Frameworks/QtDBus.framework/Versions/5/QtDBus"
	"@executable_path/../Frameworks/QtDBus.framework/Versions/5/QtDBus"
	"@executable_path/../Frameworks/QtDBus.framework/Versions/5/QtDBus"
	"@executable_path/../Frameworks/QtDBus.framework/Versions/5/QtDBus"
	"@executable_path/../Frameworks/qwt.framework/Versions/6/qwt"
	"@executable_path/../Frameworks/Python.framework/Versions/${PY_VER}/Python"
	"@executable_path/../Libraries/libMinuit2.0.dylib"
	"@executable_path/../Libraries/libboost_system-mt.dylib"
	"@executable_path/../Libraries/libboost_filesystem-mt.dylib"
	"@executable_path/../Libraries/libboost_iostreams-mt.dylib"
	"@executable_path/../Libraries/libboost_program_options-mt.dylib"
	"@executable_path/../Libraries/libboost_python39-mt.dylib"
	"@executable_path/../Libraries/libfreetype.6.dylib"
	"@executable_path/../Libraries/libpng16.16.dylib"
	"@executable_path/../Libraries/libjpeg.9.dylib"
	"@executable_path/../Libraries/libjpeg.9.dylib"
	"@executable_path/../Libraries/libtiff.5.dylib"
	"@executable_path/../Libraries/libopenblas.0.dylib"
	"@executable_path/../Libraries/libgfortran.5.dylib"
	"@executable_path/../Libraries/libgomp.1.dylib"
	"@executable_path/../Libraries/libquadmath.0.dylib"
	"@executable_path/../Libraries/libquadmath.0.dylib"
	"@executable_path/../Libraries/libgcc_s.1.dylib"
	"@executable_path/../Libraries/libgcc_s.1.dylib"
)


CNT=$(expr ${#changefrom[*]} - 1)


function fix_name()
{
	cfile=$1

	if [ ! -e ${cfile} ]; then
		echo -e "Warning: ${cfile} does not exist."
		return
	fi

	echo -e "Processing binary \"${cfile}\"..."
	chmod a+rx ${cfile}

	for idx in $(seq 0 ${CNT}); do
		cfrom=${changefrom[$idx]}
		cto=${changeto[$idx]}

		echo -e "\tChanging \"${cfrom}\"\n\t -> \"${cto}\"."
		chmod u+w ${cfile}
		${TOOL} -change ${cfrom} ${cto} ${cfile}
	done

	${STRIP} ${cfile}
	echo -e ""
}


# fix names for given files
for cfile in ${filestochange[@]}; do
	fix_name $cfile
done


# fix names for libraries
find ${PRG}/Contents/ \( -name "*.dylib" -o -name "*.so" \) -print0 \
	| while read -d $'\0' dylib; do
	fix_name $dylib
done
