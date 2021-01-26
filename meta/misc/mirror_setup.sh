#!/bin/bash
#
# automatically mirror the repositories
# @author Tobias Weber <tweber@ill.fr>
# @date 3-nov-2020
# @license GPLv2
#

mirror_from=org_repo

declare -a repos_orig=(
	"https://code.ill.fr/scientific-software/takin/meta"
	"https://code.ill.fr/scientific-software/takin/core"
	"https://code.ill.fr/scientific-software/takin/mag-core"
	"https://code.ill.fr/scientific-software/takin/tlibs"
	"https://code.ill.fr/scientific-software/takin/tlibs2"
	"https://code.ill.fr/scientific-software/takin/data"
	"https://code.ill.fr/scientific-software/takin/plugins/mnsi"
)

declare -a repos_mirror=(
	"https://github.com/tweber-ill/ill_mirror-takin2-meta"
	"https://github.com/tweber-ill/ill_mirror-takin2-core"
	"https://github.com/tweber-ill/ill_mirror-takin2-mag-core"
	"https://github.com/tweber-ill/ill_mirror-takin2-tlibs"
	"https://github.com/tweber-ill/ill_mirror-takin2-tlibs2"
	"https://github.com/tweber-ill/ill_mirror-takin2-data"
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
