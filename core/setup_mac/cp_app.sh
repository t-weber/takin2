#!/bin/bash
#
# @date 12-apr-17
# @author Tobias Weber <tobias.weber@tum.de>
# @license GPLv2
#
# copy framework libs
#
# pack with: hdiutil create takin.dmg -verbose -srcfolder takin.app -fs UDF -format "UDBZ" -volname "takin"
#

clean_frameworks=0

PRG="takin.app"

BIN_DIR="${PRG}/Contents/MacOS/"
DST_FRAMEWORK_DIR="${PRG}/Contents/Frameworks/"
DST_LIB_DIR="${PRG}/Contents/Libraries/"
DST_PLUGIN_DIR="${PRG}/Contents/PlugIns/"
DST_SITEPACKAGES_DIR="${PRG}/Contents/site-packages"
PLUGIN_DIR="/usr/local/opt/qt5/plugins/"


# -----------------------------------------------------------------------------
# frameworks
declare -a SRC_FRAMEWORKS=(
	"/usr/local/opt/qt5/lib/QtCore.framework"
	"/usr/local/opt/qt5/lib/QtGui.framework"
	"/usr/local/opt/qt5/lib/QtWidgets.framework"
	"/usr/local/opt/qt5/lib/QtOpenGL.framework"
	"/usr/local/opt/qt5/lib/QtConcurrent.framework"
	"/usr/local/opt/qt5/lib/QtXml.framework"
	"/usr/local/opt/qt5/lib/QtSvg.framework"
	"/usr/local/opt/qt5/lib/QtPrintSupport.framework"
	"/usr/local/opt/qt5/lib/QtDBus.framework"
	"/usr/local/opt/qwt/lib/qwt.framework"
	"/usr/local/opt/python/Frameworks/Python.framework"
)


# libs which need to have their symbolic link resolved
declare -a SRC_LIBS=(
	"/usr/local/opt/minuit2/lib/libMinuit2.0.dylib"
	"/usr/local/opt/boost/lib/libboost_system-mt.dylib"
	"/usr/local/opt/boost/lib/libboost_filesystem-mt.dylib"
	"/usr/local/opt/boost/lib/libboost_iostreams-mt.dylib"
	"/usr/local/opt/boost/lib/libboost_program_options-mt.dylib"
	"/usr/local/opt/boost-python3/lib/libboost_python39-mt.dylib"
	"/usr/local/opt/freetype/lib/libfreetype.6.dylib"
	"/usr/local/opt/libpng/lib/libpng16.16.dylib"
	"/usr/local/opt/libjpeg/lib/libjpeg.9.dylib"
	"/usr/local/opt/libtiff/lib/libtiff.5.dylib"
	"/usr/local/opt/gcc/lib/gcc/10/libgfortran.5.dylib"
	"/usr/local/opt/gcc/lib/gcc/10/libquadmath.0.dylib"
	"/usr/local/opt/gcc/lib/gcc/10/libgomp.1.dylib"
	"/usr/local/opt/gcc/lib/gcc/10/libgcc_s.1.dylib"
	"/usr/local/opt/openblas/lib/libopenblas.0.dylib"
)


# qt plugins
declare -a SRC_PLUGINS=(
	"printsupport/libcocoaprintersupport.dylib"
	"imageformats/libqsvg.dylib"
	"imageformats/libqicns.dylib"
	"imageformats/libqjpeg.dylib"
	"iconengines/libqsvgicon.dylib"
	"styles/libqmacstyle.dylib"
	"platforms/libqcocoa.dylib"
	"platforms/libqminimal.dylib"
)
# -----------------------------------------------------------------------------



# -----------------------------------------------------------------------------
# clean old files
rm -rfv ${BIN_DIR}
rm -rfv ${DST_FRAMEWORK_DIR}
rm -rfv ${DST_LIB_DIR}
rm -rfv ${DST_PLUGIN_DIR}
rm -rfv ${DST_SITEPACKAGES_DIR}
# -----------------------------------------------------------------------------



# -----------------------------------------------------------------------------
# create dirs
mkdir -pv "${DST_FRAMEWORK_DIR}"
mkdir -pv "${DST_LIB_DIR}"
mkdir -pv "${DST_PLUGIN_DIR}/printsupport"
mkdir -pv "${DST_PLUGIN_DIR}/imageformats"
mkdir -pv "${DST_PLUGIN_DIR}/iconengines"
mkdir -pv "${DST_PLUGIN_DIR}/platforms"
mkdir -pv "${DST_PLUGIN_DIR}/platformthemes"
mkdir -pv "${DST_PLUGIN_DIR}/styles"
mkdir -pv "${BIN_DIR}"
mkdir -pv "${PRG}/Contents/Resources"
mkdir -pv "${DST_SITEPACKAGES_DIR}"
# -----------------------------------------------------------------------------



# -----------------------------------------------------------------------------
# copy frameworks
for srclib in ${SRC_FRAMEWORKS[@]}; do
	echo -e "Copying framework \"${srclib}\"..."
	cp -rv ${srclib} "${DST_FRAMEWORK_DIR}"
done


# copy libs, following links
for srclib in ${SRC_LIBS[@]}; do
	echo -e "Copying library \"${srclib}\"..."
	cp -Lrv ${srclib} "${DST_LIB_DIR}"
done


# copy plugins
for srclib in ${SRC_PLUGINS[@]}; do
	echo -e "Copying plugin \"${srclib}\"..."
	cp -v "${PLUGIN_DIR}${srclib}" "${DST_PLUGIN_DIR}${srclib}"
	chmod a+rx "${DST_PLUGIN_DIR}${srclib}"
done


# copy binaries
cp -v bin/takin "${BIN_DIR}"
cp -v bin/takin_cif2xml "${BIN_DIR}"
cp -v bin/takin_findsg "${BIN_DIR}"
cp -v bin/takin_pol "${BIN_DIR}"
cp -v bin/takinmod_py "${BIN_DIR}"
cp -v bin/takinmod_jl "${BIN_DIR}"
cp -v bin/takin_structfact "${BIN_DIR}"
cp -v bin/takin_magstructfact "${BIN_DIR}"
cp -v bin/takin_scanbrowser "${BIN_DIR}"
cp -v bin/takin_magsgbrowser "${BIN_DIR}"
cp -v bin/takin_moldyn "${BIN_DIR}"

cp -v bin/takin_convofit "${BIN_DIR}"
cp -v bin/takin_convoseries "${BIN_DIR}"
cp -v bin/takin_polextract "${BIN_DIR}"


# data files
cp -v setup_mac/Info.plist "${PRG}/Contents/"
cp -rv data/res "${PRG}/Contents/"
cp -rv doc/* "${PRG}/Contents/res/doc/"
cp -v data/res/icons/takin.icns "${PRG}/Contents/Resources/"
cp -v *.txt "${PRG}/Contents/Resources/"
cp -rv 3rdparty_licenses "${PRG}/Contents/Resources/"
rm -v "${PRG}/Contents/Resources/CMakeLists.txt"


# python site-packages
cp -rv /usr/local/opt/numpy/lib/python3.9/site-packages/numpy "${DST_SITEPACKAGES_DIR}"
cp -rv /usr/local/opt/scipy/lib/python3.9/site-packages/scipy "${DST_SITEPACKAGES_DIR}"
#ln -sf "../site-packages" "${BIN_DIR}/site-packages"
#ln -sf ./Frameworks/Python.framework/Versions/Current/lib/python3.9/site-packages ${DST_SITEPACKAGES_DIR}
# -----------------------------------------------------------------------------



# -----------------------------------------------------------------------------
# attributes
chmod -R u+rw,a+r "${DST_FRAMEWORK_DIR}"
chmod -R u+rw,a+r "${DST_LIB_DIR}"
chmod -R u+rw,a+r "${DST_PLUGIN_DIR}"

# only the binaries in BIN_DIR should have the executable flag
find "${PRG}" -type f -print0 | xargs -0 chmod a-x
chmod a+x ${BIN_DIR}*

find "${DST_FRAMEWORK_DIR}" -type d -print0 | xargs -0 chmod a+x
find "${DST_LIB_DIR}" -type d -print0 | xargs -0 chmod a+x
find "${DST_PLUGIN_DIR}" -type d -print0 | xargs -0 chmod a+x
# -----------------------------------------------------------------------------



# -----------------------------------------------------------------------------
# delete unnecessary files
if [ $clean_frameworks -ne 0 ]; then
	# clean site-packages
	rm -rfv ${DST_SITEPACKAGES_DIR}/PyQt5*
	rm -rfv ${DST_SITEPACKAGES_DIR}/wheel*
	rm -rfv ${DST_SITEPACKAGES_DIR}/setuptools*
	rm -rfv ${DST_SITEPACKAGES_DIR}/pip*
	rm -rfv ${DST_SITEPACKAGES_DIR}/distutils*
	rm -rfv ${DST_SITEPACKAGES_DIR}/_distutils*
	rm -rfv ${DST_SITEPACKAGES_DIR}/__pycache*
	rm -rfv ${DST_SITEPACKAGES_DIR}/pkg_*
	rm -rfv ${DST_SITEPACKAGES_DIR}/sip*
	rm -rfv ${DST_SITEPACKAGES_DIR}/site*
	rm -rfv ${DST_SITEPACKAGES_DIR}/easy_*

	# clean non-needed files from frameworks
	find "${DST_FRAMEWORK_DIR}" -type d -name "Headers" -print0 | xargs -0 rm -rfv
	find "${DST_FRAMEWORK_DIR}" -type d -name "Current" -print0 | xargs -0 rm -rfv

	declare -a QT_FW_LIBS=("QtCore" "QtGui" "QtWidgets" "QtConcurrent" "QtOpenGL"
		"QtSvg" "QtXml" "QtDBus" "QtPrintSupport" "qwt")

	for qlib in ${QT_FW_LIBS[@]}; do
		rm -rfv ${DST_FRAMEWORK_DIR}/${qlib}.framework/Resources
		rm -rfv ${DST_FRAMEWORK_DIR}/${qlib}.framework/${qlib}*
	done
fi
# -----------------------------------------------------------------------------
