#!/bin/bash
#
# cleans temporary files
# @author Tobias Weber <tobias.weber@tum.de>
# @license GPLv2
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

find bin -regex 'bin/[_a-zA-Z0-9]*' | xargs rm -fv
rm -fv bin/*.exe
rm -fv plugins/*.so
rm -fv plugins/*.dll

rm -fv doc/takin.qch
rm -fv doc/takin.qhc
rm -rfv doc/devel/html
rm -fv doc/devel/*.tmp

if [ -f Makefile ]
then
	echo -e "\nCleaning stuff made by Makefile..."
	make clean
fi

rm -fv CMakeCache.txt
rm -rfv CMakeFiles

rm -rfv build/
