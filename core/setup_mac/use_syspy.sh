#!/bin/bash
#
# @date 16-dec-2020
# @author Tobias Weber <tweber@ill.fr>
# @license GPLv2
#
# links to the system python libraries instead and removes the packaged ones
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

TOOL=install_name_tool
STRIP=strip

PY_VER=3.9


# files whose linkage is to be changed
declare -a filestochange=(
	"${PRG}/Contents/MacOS/takinmod_py"
)



# original symbols
declare -a changefrom=(
	"@executable_path/../Frameworks/Python.framework/Versions/${PY_VER}/Python"
)


# symbols to change into
declare -a changeto=(
	"/Library/Frameworks/Python.framework/Versions/Current/Python"
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


# remove the packaged python frameworks and now non-needed libraries
rm -rfv "${PRG}/Contents/Frameworks/Python.framework"
rm -rfv "${PRG}/Contents/site-packages"
rm -fv "${PRG}/Contents/Libraries/libgfortran.5.dylib"
rm -fv "${PRG}/Contents/Libraries/libquadmath.0.dylib"
rm -fv "${PRG}/Contents/Libraries/libgomp.1.dylib"
rm -fv "${PRG}/Contents/Libraries/libgcc_s.1.dylib"
rm -fv "${PRG}/Contents/Libraries/libopenblas.0.dylib"
