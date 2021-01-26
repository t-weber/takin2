#!/bin/bash
#
# @date 8-dec-20
# @author Tobias Weber <tweber@ill.fr>
# @license GPLv2
#
# finds binaries with bindings to libraries in /usr/local
#

PRG="takin.app"


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


function find_bindings()
{
	cfile=$1

	if [ ! -e ${cfile} ]; then
		echo -e "Warning: ${cfile} does not exist."
		return
	fi

	echo -e "--------------------------------------------------------------------------------"
	echo -e "\"${cfile}\":"
	otool -L ${cfile} | grep --color /usr/local
	echo -e "--------------------------------------------------------------------------------\n"
}


# find local bindings for given files
for cfile in ${filestochange[@]}; do
	find_bindings $cfile
done


# find local bindings for libraries
find ${PRG}/Contents/ \( -name "*.dylib" -o -name "*.so" \) -print0 \
	| while read -d $'\0' dylib; do
	find_bindings $dylib
done
