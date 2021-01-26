#!/bin/bash
#
# cleans temporary files
# @author Tobias Weber <tobias.weber@tum.de>
# @license GPLv2
#

find bin -regex 'bin/[_a-zA-Z0-9]*' | xargs rm -fv
rm -fv bin/*.exe
rm -fv plugins/*.so
rm -fv plugins/*.dll

rm -fv doc/takin.qch
rm -fv doc/takin.qhc
rm -rfv doc/devel/html
rm -fv doc/devel/*.tmp

if [ -f Makefile ]
then
	echo -e "\nCleaning stuff made by Makefile..."
	make clean
fi

rm -fv CMakeCache.txt
rm -rfv CMakeFiles

rm -rfv build/
