#!/bin/bash
#
# install all packages needed for building
# @author Tobias Weber <tobias.weber@tum.de>
# @date 21-jul-22
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

#if [[ $(id -u) -gt 0 ]]; then
#	echo -e "Please run this script as root."
#	exit -1
#fi


# -----------------------------------------------------------------------------
# install packages
# -----------------------------------------------------------------------------
if ! sudo pacman -S cmake gcc boost \
	qt5 qt5-webkit qwt qcustomplot \
	qhull lapacke gnuplot python3
then
	echo -e "Error: Could not install packages necessary for building."
	exit -1
fi
# -----------------------------------------------------------------------------
