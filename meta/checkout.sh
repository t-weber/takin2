#!/bin/bash
#
# checkout all of takin's repos
# @author Tobias Weber <tweber@ill.fr>
# @date jul-20
# @license GPLv2
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

${GIT} clone https://code.ill.fr/scientific-software/takin/meta.git
${GIT} clone https://code.ill.fr/scientific-software/takin/core.git
${GIT} clone https://code.ill.fr/scientific-software/takin/mag-core.git
${GIT} clone https://code.ill.fr/scientific-software/takin/tlibs.git
${GIT} clone https://code.ill.fr/scientific-software/takin/tlibs2.git
${GIT} clone https://code.ill.fr/scientific-software/takin/data.git

mkdir plugins
pushd plugins
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
popd



echo -e "\n--------------------------------------------------------------------------------"
echo -e "Done, the Takin source can be found in \"${TAKIN_DIR}\"."
echo -e "--------------------------------------------------------------------------------\n"
