#!/bin/bash
#
# automatically mirror the repositories
# @author Tobias Weber <tweber@ill.fr>
# @date sep-2020
# @license GPLv2


mirror_from=org_repo


declare -a repo_dirs=(
	"ill_mirror-takin2-meta"
	"ill_mirror-takin2-core"
	"ill_mirror-takin2-mag-core"
	"ill_mirror-takin2-tlibs"
	"ill_mirror-takin2-tlibs2"
	"ill_mirror-takin2-data"
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
