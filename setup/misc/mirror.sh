#!/bin/bash
#
# automatically mirror the repositories
# @author Tobias Weber <tweber@ill.fr>
# @date sep-2020
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

mirror_from=org_repo


declare -a repo_dirs=(
	"ill_mirror-takin2-meta"
	"ill_mirror-takin2-core"
	"ill_mirror-takin2-mag-core"
	"ill_mirror-takin2-tlibs"
	"ill_mirror-takin2-tlibs2"
	"ill_mirror-takin2-data"
	"ill_mirror-takin2-paths"
	"ill_mirror-takin2-mnsi"
)


for repo_dir in ${repo_dirs[@]}; do
	echo -e "--------------------------------------------------------------------------------"
	echo -e "Mirroring ${repo_dir}"
	echo -e "--------------------------------------------------------------------------------"

	pushd "${repo_dir}"
	git pull -v "${mirror_from}" master
	git push -v
	popd

	echo -e "--------------------------------------------------------------------------------\n"
done
