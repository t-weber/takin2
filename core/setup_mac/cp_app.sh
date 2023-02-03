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
# ----------------------------------------------------------------------------
# Takin (inelastic neutron scattering software package)
# Copyright (C) 2017-2023  Tobias WEBER (Institut Laue-Langevin (ILL),
#                          Grenoble, France).
# Copyright (C) 2013-2017  Tobias WEBER (Technische Universitaet Muenchen
#                          (TUM), Garching, Germany).
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; version 2 of the License.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along
# with this program; if not, write to the Free Software Foundation, Inc.,
# 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
# ----------------------------------------------------------------------------
#

clean_frameworks=1

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
	"/usr/local/opt/qwt-qt5/lib/qwt.framework"
	"/Library/Frameworks/Python.framework"
)
#	"/usr/local/opt/python/Frameworks/Python.framework"

# libs which need to have their symbolic link resolved
declare -a SRC_LIBS=(
	"/usr/local/opt/minuit2/lib/libMinuit2.0.dylib"
	"/usr/local/opt/boost/lib/libboost_system-mt.dylib"
	"/usr/local/opt/boost/lib/libboost_filesystem-mt.dylib"
	"/usr/local/opt/boost/lib/libboost_atomic-mt.dylib"
	"/usr/local/opt/boost/lib/libboost_iostreams-mt.dylib"
	"/usr/local/opt/boost/lib/libboost_program_options-mt.dylib"
	"/usr/local/opt/boost-python3/lib/libboost_python311-mt.dylib"
	"/usr/local/opt/freetype/lib/libfreetype.6.dylib"
	"/usr/local/opt/libpng/lib/libpng16.16.dylib"
	"/usr/local/opt/libjpeg/lib/libjpeg.9.dylib"
	"/usr/local/opt/libtiff/lib/libtiff.5.dylib"
	"/usr/local/opt/lapack/lib/liblapacke.3.dylib"
	"/usr/local/opt/lapack/lib/liblapack.3.dylib"
	"/usr/local/opt/lapack/lib/libblas.3.dylib"
	"/usr/local/opt/gcc/lib/gcc/current/libgfortran.5.dylib"
	"/usr/local/opt/gcc/lib/gcc/current/libquadmath.0.dylib"
	"/usr/local/opt/gcc/lib/gcc/current/libgomp.1.dylib"
	"/usr/local/opt/gcc/lib/gcc/current/libgcc_s.1.1.dylib"
	"/usr/local/opt/glib/lib/libgthread-2.0.0.dylib"
	"/usr/local/opt/glib/lib/libglib-2.0.0.dylib"
	"/usr/local/opt/pcre/lib/libpcre.1.dylib"
	"/usr/local/opt/pcre2/lib/libpcre2-16.0.dylib"
	"/usr/local/opt/pcre2/lib/libpcre2-8.0.dylib"
	"/usr/local/opt/gettext/lib/libintl.8.dylib"
	"/usr/local/opt/zstd/lib/libzstd.1.dylib"
	"/usr/local/opt/xz/lib/liblzma.5.dylib"
	"/usr/local/opt/qhull/lib/libqhull_r.8.0.dylib"
	"/usr/local/opt/hdf5/lib/libhdf5_cpp.200.dylib"
	"/usr/local/opt/hdf5/lib/libhdf5.200.dylib"
	"/usr/local/opt/libaec/lib/libsz.2.dylib"
	"/usr/local/lib/libqcustomplot.dylib"
)
#	"/usr/local/opt/openblas/lib/libopenblas.0.dylib"


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
cp -v bin/takin_bz "${BIN_DIR}"
cp -v bin/takinmod_py "${BIN_DIR}"
cp -v bin/takinmod_jl "${BIN_DIR}"
cp -v bin/takin_structfact "${BIN_DIR}"
cp -v bin/takin_magstructfact "${BIN_DIR}"
cp -v bin/takin_scanbrowser "${BIN_DIR}"
cp -v bin/takin_magsgbrowser "${BIN_DIR}"
cp -v bin/takin_magdyn "${BIN_DIR}"
cp -v bin/takin_moldyn "${BIN_DIR}"

cp -v bin/takin_convofit "${BIN_DIR}"
cp -v bin/takin_convoseries "${BIN_DIR}"
cp -v bin/takin_polextract "${BIN_DIR}"

# copy plugin modules
cp -v plugins/*.dylib "${DST_PLUGIN_DIR}"


# data files
cp -v setup_mac/Info.plist "${PRG}/Contents/"
cp -rv data/res "${PRG}/Contents/"
cp -rv doc/* "${PRG}/Contents/res/doc/"
cp -v data/res/icons/takin.icns "${PRG}/Contents/Resources/"
cp -v *.txt "${PRG}/Contents/Resources/"
cp -rv 3rdparty_licenses "${PRG}/Contents/Resources/"
rm -v "${PRG}/Contents/Resources/CMakeLists.txt"


# python site-packages
#cp -rv /usr/local/opt/numpy/lib/python3.9/site-packages/numpy "${DST_SITEPACKAGES_DIR}"
#cp -rv /usr/local/opt/scipy/lib/python3.9/site-packages/scipy "${DST_SITEPACKAGES_DIR}"
##ln -sf "../site-packages" "${BIN_DIR}/site-packages"
##ln -sf ./Frameworks/Python.framework/Versions/Current/lib/python3.9/site-packages ${DST_SITEPACKAGES_DIR}
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
	rm -rfv ${DST_SITEPACKAGES_DIR}/__pycache__*
	rm -rfv ${DST_SITEPACKAGES_DIR}/pkg_*
	rm -rfv ${DST_SITEPACKAGES_DIR}/sip*
	rm -rfv ${DST_SITEPACKAGES_DIR}/site*
	rm -rfv ${DST_SITEPACKAGES_DIR}/easy_*

	# clean non-needed files from frameworks
	find "${DST_FRAMEWORK_DIR}" -type d -name "Headers" -print0 | xargs -0 rm -rfv
	find "${DST_FRAMEWORK_DIR}" -type d -name "Current" -print0 | xargs -0 rm -rfv
	find "${DST_FRAMEWORK_DIR}" -name "*.pyc" -print0 | xargs -0 rm -rfv

	declare -a QT_FW_LIBS=("QtCore" "QtGui" "QtWidgets" "QtConcurrent" "QtOpenGL"
		"QtSvg" "QtXml" "QtDBus" "QtPrintSupport" "qwt")

	for qlib in ${QT_FW_LIBS[@]}; do
		rm -rfv ${DST_FRAMEWORK_DIR}/${qlib}.framework/Resources
		rm -rfv ${DST_FRAMEWORK_DIR}/${qlib}.framework/Headers
		rm -rfv ${DST_FRAMEWORK_DIR}/${qlib}.framework/${qlib}*
	done

	# remove old signatures
	codesign --remove-signature ${DST_FRAMEWORK_DIR}/Python.framework/Versions/Current/Python
fi


find ${PRG} -type f -name ".dir" -exec rm -fv {} \; -print
find ${PRG} -type f -name ".DS_Store" -exec rm -fv {} \; -print
# -----------------------------------------------------------------------------
