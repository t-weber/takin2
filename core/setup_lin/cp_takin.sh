#!/bin/bash
#
# creates a distro archive
# @author Tobias Weber <tobias.weber@tum.de>
# @date 2016 - 2023
# @license GPLv2
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


# installation directory
INSTDIR="$1"

if [ "${INSTDIR}" = "" ]; then
	INSTDIR=takin
fi

mkdir -p ${INSTDIR}


# main programs
cp -rv bin			${INSTDIR}/
cp -v takin.sh			${INSTDIR}/


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


# remove uneeded files
find ${INSTDIR}/ -type f -name ".dir" -exec rm -v {} \; -print

# stripping
strip -v ${INSTDIR}/bin/*
