#!/bin/bash
#
# @date 20-jan-21
# @author Tobias Weber <tweber@ill.fr>
# @license GPLv2
#
# finds library dependencies
#

if [ $# -lt 2 ]; then
	echo -e "Please give a path and a symbol."
	exit -1
fi

dir=$1
symtofind=$2

find ${dir} \( -name "*.dylib" -o -name "*.so" \) -print0 \
	| while read -d $'\0' dylib; do

	if [ ! -e ${dylib} ]; then
		echo -e "Warning: \"${dylib}\" does not exist."
		continue
	fi

	result=$(otool -L ${dylib} | grep -i --color ${symtofind})
	if [ "${result}" != "" ]; then
		echo -e "Library: ${dylib}:\n${result}\n"
	fi
done
