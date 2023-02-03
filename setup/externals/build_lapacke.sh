#!/bin/bash
#
# builds lapack and lapacke
# @author Tobias Weber <tweber@ill.fr>
# @date jul-2022
# @license GPLv2
#
# ----------------------------------------------------------------------------
# Takin (inelastic neutron scattering software package)
# Copyright (C) 2017-2022  Tobias WEBER (Institut Laue-Langevin (ILL),
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

NUM_CORES=$(nproc)


BUILD_FOR_MINGW=0
if [ "$1" == "--mingw" ]; then
	BUILD_FOR_MINGW=1
fi



LAPACK_REMOTE=https://codeload.github.com/Reference-LAPACK/lapack/zip/refs/heads/master
LAPACK_LOCAL=${LAPACK_REMOTE##*[/\\]}


rm -f "${LAPACK_LOCAL}"


if ! wget ${LAPACK_REMOTE}; then
	echo -e "Could not download ${LAPACK_REMOTE}."
	exit -1
fi


unzip "${LAPACK_LOCAL}"


cd lapack-master
mkdir build && cd build


if [ $BUILD_FOR_MINGW -ne 0 ]; then
	mingw64-cmake -DCMAKE_BUILD_TYPE=Release -DLAPACKE=TRUE ..
	mingw64-make -j${NUM_CORES} && sudo mingw64-make install
else
	cmake -DCMAKE_BUILD_TYPE=Release -DLAPACKE=TRUE -DBUILD_SHARED_LIBS=TRUE ..
	make -j${NUM_CORES} && sudo make install
fi
