#!/bin/bash
# prebuilding steps
# @author Tobias Weber <tobias.weber@tum.de>
# @license GPLv2
#


if [[ $# -ge 2 ]]; then
	COLGEN="$(dirname $2)/qhelpgenerator"
	COLGEN_ALT="$(dirname $2)/qcollectiongenerator"

	TAKINROOT="$1"
else
	COLGEN="$(which qhelpgenerator-qt5)"
	if [ $? -ne 0 ]; then
		COLGEN="$(which qhelpgenerator)"
	fi

	COLGEN_ALT="$(which qcollectiongenerator-qt5)"
	if [ $? -ne 0 ]; then
		COLGEN_ALT="$(which qcollectiongenerator)"
	fi

	TAKINROOT=.
fi


echo -e "Takin root dir: ${TAKINROOT}\n"



echo -e "--------------------------------------------------------------------------------"
echo -e "building docs using collection generator: ${COLGEN}..."

${COLGEN} ${TAKINROOT}/doc/takin.qhcp -o ${TAKINROOT}/doc/takin.qhc

# try alternate help collection generator
if [ ! -f ${TAKINROOT}/doc/takin.qhc ]; then
	${COLGEN_ALT} ${TAKINROOT}/doc/takin.qhcp -o ${TAKINROOT}/doc/takin.qhc
fi

cp -v ${TAKINROOT}/doc/takin.qhc ${TAKINROOT}/res/doc/
cp -v ${TAKINROOT}/doc/takin.qch ${TAKINROOT}/res/doc/

echo -e "--------------------------------------------------------------------------------\n"



echo -e "--------------------------------------------------------------------------------"
echo -e "building devel docs..."

cd ${TAKINROOT}
doxygen takin-doc

echo -e "--------------------------------------------------------------------------------\n"

