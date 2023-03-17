#!/bin/bash
#
# automatically mirror the repositories
# @author Tobias Weber <tweber@ill.fr>
# @date 3-nov-2020
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

declare -a repos_orig=(
	"https://code.ill.fr/scientific-software/takin/setup"
	"https://code.ill.fr/scientific-software/takin/core"
	"https://code.ill.fr/scientific-software/takin/mag-core"
	"https://code.ill.fr/scientific-software/takin/tlibs"
	"https://code.ill.fr/scientific-software/takin/tlibs2"
	"https://code.ill.fr/scientific-software/takin/data"
	"https://code.ill.fr/scientific-software/takin/paths"
	"https://code.ill.fr/scientific-software/takin/plugins/mnsi"
)

declare -a repos_mirror=(
	"https://github.com/tweber-ill/ill_mirror-takin2-setup"
	"https://github.com/tweber-ill/ill_mirror-takin2-core"
	"https://github.com/tweber-ill/ill_mirror-takin2-mag-core"
	"https://github.com/tweber-ill/ill_mirror-takin2-tlibs"
	"https://github.com/tweber-ill/ill_mirror-takin2-tlibs2"
	"https://github.com/tweber-ill/ill_mirror-takin2-data"
	"https://github.com/tweber-ill/ill_mirror-takin2-paths"
	"https://github.com/tweber-ill/ill_mirror-takin2-mnsi"
)


mkdir -p mirror
cd mirror

for (( i=0; i<${#repos_orig[@]}; ++i )); do
    repo_orig=${repos_orig[$i]}
    repo_mirror=${repos_mirror[$i]}
    repo_dir=${repo_mirror##*[/\\]}

	echo -e "--------------------------------------------------------------------------------"
	echo -e "Mirroring:     ${repo_dir}"
	echo -e "Original repo: ${repo_orig}"
	echo -e "Mirror repo:   ${repo_mirror}"
	echo -e "--------------------------------------------------------------------------------"

    git clone ${repo_mirror}

    pushd "${repo_dir}"
    git remote add ${mirror_from} ${repo_orig}
	popd

	echo -e "--------------------------------------------------------------------------------\n"
done
