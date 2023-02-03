#!/bin/bash
#
# @date 20-jan-21
# @author Tobias Weber <tweber@ill.fr>
# @license GPLv2
#
# finds library dependencies
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

if [ $# -lt 2 ]; then
	echo -e "Please give a path and a symbol."
	exit -1
fi

dir=$1
symtofind=$2

find ${dir} \( -name "*.dylib" -o -name "*.so" \) -print0 \
	| while read -d $'\0' dylib; do

	if [ ! -e ${dylib} ]; then
		echo -e "Warning: \"${dylib}\" does not exist."
		continue
	fi

	result=$(otool -L ${dylib} | grep -i --color ${symtofind})
	if [ "${result}" != "" ]; then
		echo -e "Library: ${dylib}:\n${result}\n"
	fi
done
