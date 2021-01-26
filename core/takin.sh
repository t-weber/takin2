#!/bin/sh
#
# calling takin with local libs
# @author Tobias Weber <tobias.weber@tum.de>
# @license GPLv2
#

TAKINDIR=$(dirname $0)
echo "Takin directory: ${TAKINDIR}"


# look for binary in bin folder
TAKINBIN=${TAKINDIR}/bin/takin

if [ ! -f ${TAKINBIN} ]; then
	# look for binary in root folder
	TAKINBIN=${TAKINDIR}/takin

	if [ ! -f ${TAKINBIN} ]; then
		echo -e "Takin binary was not found!"
		exit -1
	fi
fi

LD_LIBRARY_PATH=./lib:${TAKINDIR}/lib:/usr/local/lib:/usr/local/lib64:$LD_LIBRARY_PATH \
	${TAKINBIN} "$@"

