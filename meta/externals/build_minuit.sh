#!/bin/bash
#
# builds minuit
# @author Tobias Weber <tweber@ill.fr>
# @date sep-2020
# @license GPLv2
#

NUM_CORES=$(nproc)


BUILD_FOR_MINGW=0
if [ "$1" == "--mingw" ]; then
	BUILD_FOR_MINGW=1
fi



MINUIT_REMOTE=http://www.cern.ch/mathlibs/sw/5_34_14/Minuit2/Minuit2-5.34.14.tar.gz
MINUIT_LOCAL=${MINUIT_REMOTE##*[/\\]}
MINUIT_DIR=${MINUIT_LOCAL%.tar.gz}


rm -f "${MINUIT_LOCAL}"
rm -rf "${MINUIT_DIR}"


if ! wget ${MINUIT_REMOTE}; then
	echo -e "Could not download ${MINUIT_REMOTE}."
	exit -1
fi


tar xvf "${MINUIT_LOCAL}"


cd "${MINUIT_DIR}"


if BUILD_FOR_MINGW; then
	mingw64-configure --disable-openmp
	mingw64-make -j${NUM_CORES} && sudo mingw64-make install-strip
else
	./configure --disable-openmp
	make -j${NUM_CORES} && sudo make install-strip
fi
