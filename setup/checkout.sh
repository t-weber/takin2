#!/bin/bash
#
# checkout all of takin's repos
# @author Tobias Weber <tweber@ill.fr>
# @date jul-20
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

echo -e "================================================================================"
echo -e "Setting up Takin 2 sources"
echo -e "================================================================================"



echo -e "\n--------------------------------------------------------------------------------"
echo -e "Preliminary steps..."
echo -e "--------------------------------------------------------------------------------"

GIT="$(which git)"
if [ "$GIT" = "" ]; then
	echo -e "\nERROR: The \"git\" tool was not found, please install it before running this script.\n"
	exit -1
fi


TAKIN_DIR=takin
if [ -L ${TAKIN_DIR} ] || [ -d ${TAKIN_DIR} ]; then
	echo -e "\nERROR: A \"${TAKIN_DIR}\" directory already exists here, exiting.\n"
	exit -1
fi

mkdir -v ${TAKIN_DIR}
cd ${TAKIN_DIR}



echo -e "\n--------------------------------------------------------------------------------"
echo -e "Cloning all Takin repositories..."
echo -e "--------------------------------------------------------------------------------"

${GIT} clone https://code.ill.fr/scientific-software/takin/setup.git
${GIT} clone https://code.ill.fr/scientific-software/takin/core.git
${GIT} clone https://code.ill.fr/scientific-software/takin/mag-core.git
${GIT} clone https://code.ill.fr/scientific-software/takin/tlibs.git
${GIT} clone https://code.ill.fr/scientific-software/takin/tlibs2.git
${GIT} clone https://code.ill.fr/scientific-software/takin/data.git
${GIT} clone https://code.ill.fr/scientific-software/takin/paths.git

mkdir plugins
pushd plugins
${GIT} clone https://code.ill.fr/scientific-software/takin/plugins/magnons.git
${GIT} clone https://code.ill.fr/scientific-software/takin/plugins/mnsi.git
popd



echo -e "\n--------------------------------------------------------------------------------"
echo -e "Creating links..."
echo -e "--------------------------------------------------------------------------------"

pushd core
ln -sf ../data
ln -sf ../tlibs
popd

pushd mag-core
ln -sf ../data
ln -sf ../tlibs2
ln -sf ../paths/src pathslib
popd

pushd paths
ln -sf ../tlibs2
popd

pushd plugins
pushd magnons
pushd ext
ln -sf ../../../tlibs2
ln -sf ../../../tlibs
ln -sf ../../../core takin
ln -sf ../../../mag-core takin2
popd
popd
popd


echo -e "\n--------------------------------------------------------------------------------"
echo -e "Done, the Takin source can be found in \"${TAKIN_DIR}\"."
echo -e "--------------------------------------------------------------------------------\n"
