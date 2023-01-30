#!/bin/bash
#
# @date 16-dec-2020
# @author Tobias Weber <tweber@ill.fr>
# @license GPLv2
#
# cleans up the python distribution
#
# ----------------------------------------------------------------------------
# Takin (inelastic neutron scattering software package)
# Copyright (C) 2017-2021  Tobias WEBER (Institut Laue-Langevin (ILL),
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

PRG="takin.app"
aggressive_cleaning=1


# remove cache, object and hidden files
find ${PRG}/Contents/Frameworks/Python.framework -type d -name "__pycache__" -exec rm -rfv {} \;
find ${PRG}/Contents/Frameworks/Python.framework -type f -name "*.pyc" -exec rm -fv {} \;
find ${PRG}/Contents/Frameworks/Python.framework -type f -name "*.o" -exec rm -fv {} \;
find ${PRG}/Contents/Frameworks/Python.framework -type f -name ".*" -exec rm -fv {} \;


# remove tcl/tk
rm -fv ${PRG}/Contents/Frameworks/Python.framework/Versions/Current/lib/libtk*
rm -fv ${PRG}/Contents/Frameworks/Python.framework/Versions/Current/lib/libtcl*
rm -rfv ${PRG}/Contents/Frameworks/Python.framework/Versions/Current/lib/tcl8*
rm -rfv ${PRG}/Contents/Frameworks/Python.framework/Versions/Current/lib/tk8*
rm -rfv ${PRG}/Contents/Frameworks/Python.framework/Versions/Current/lib/itcl4*
rm -fv ${PRG}/Contents/Frameworks/Python.framework/Versions/Current/lib/pkgconfig/tcl.pc
rm -fv ${PRG}/Contents/Frameworks/Python.framework/Versions/Current/lib/tcl*
rm -fv ${PRG}/Contents/Frameworks/Python.framework/Versions/Current/lib/tcl*
rm -fv ${PRG}/Contents/Frameworks/Python.framework/Versions/Current/lib/tk*
rm -fv ${PRG}/Contents/Frameworks/Python.framework/Versions/Current/lib/pkgconfig/tk.pc
rm -rfv ${PRG}/Contents/Frameworks/Python.framework/Versions/Current/lib/python3.11/tkinter
rm -fv ${PRG}/Contents/Frameworks/Python.framework/Versions/Current/lib/python3.11/lib-dynload/_tkinter*.so
rm -fv ${PRG}/Contents/Frameworks/Python.framework/Versions/Current/lib/Tk*


# remove non-needed site packages
pushd ${PRG}/Contents
ln -sf Frameworks/Python.framework/Versions/Current/lib/python3.11/site-packages
popd

rm -rfv ${PRG}/Contents/site-packages/setuptools*
rm -rfv ${PRG}/Contents/site-packages/pip*
rm -rfv ${PRG}/Contents/site-packages/easy*
rm -rfv ${PRG}/Contents/site-packages/pkg_resources


# strip libraries
find ${PRG}/Contents/Frameworks/Python.framework -name "*.so" -exec strip -v {} \; -print
find ${PRG}/Contents/Frameworks/Python.framework -name "*.dylib" -exec strip -v {} \; -print


# remove app
rm -rfv ${PRG}/Contents/Frameworks/Python.framework/Versions/Current/Resources/Python.app

# remove tests
rm -rfv ${PRG}/Contents/Frameworks/Python.framework/Versions/Current/lib/python3.11/test


if [ $aggressive_cleaning -ne 0 ]; then
	# remove everything except Python and lib
	rm -rfv ${PRG}/Contents/Frameworks/Python.framework/Versions/Current/bin
	rm -rfv ${PRG}/Contents/Frameworks/Python.framework/Versions/Current/etc
	rm -rfv ${PRG}/Contents/Frameworks/Python.framework/Versions/Current/share
	rm -rfv ${PRG}/Contents/Frameworks/Python.framework/Versions/Current/include
	rm -rfv ${PRG}/Contents/Frameworks/Python.framework/Versions/Current/Headers
	rm -rfv ${PRG}/Contents/Frameworks/Python.framework/Versions/Current/_CodeSignature
	rm -rfv ${PRG}/Contents/Frameworks/Python.framework/Versions/Current/Resources
	rm -fv ${PRG}/Contents/Frameworks/Python.framework/Resources
	rm -fv ${PRG}/Contents/Frameworks/Python.framework/Headers

	for file in ${PRG}/Contents/Frameworks/Python.framework/Versions/Current/lib/*; do
		# remove everything except file libpython3.11.dylib and dir python3.11
		if [[ $file =~ [a-zA-Z0-9\./]*python[a-zA-Z0-9\./]* ]]; then
			echo -e "Keeping \"${file}\"."
			continue
		fi

		rm -rfv $file
	done
fi
