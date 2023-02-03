#!/bin/bash
#
# downloads external dependencies
# @author Tobias Weber <tweber@ill.fr>
# @date 9-apr-2020
# @license GPLv2, see 'LICENSE' file
#

use_archived_copies=0


if [ -L takin ] || [ -d takin ]; then
	echo -e "A takin directory already exists. Skipping.";
else
	if [ -d archived/takin ] && [ $use_archived_copies -eq 1 ]; then
		echo -e "Using local archived copy of takin.";
		ln -sf archived/takin
	else
		echo -e "Cloning takin/core..."
		git clone https://code.ill.fr/scientific-software/takin/core.git
		mv -v core takin
	fi
fi


if [ -L tlibs ] || [ -d tlibs ]; then
	echo -e "A tlibs directory already exists. Skipping.";
else
	if [ -d archived/tlibs ] && [ $use_archived_copies -eq 1 ]; then
		echo -e "Using local archived copy of tlibs.";
		ln -sf archived/tlibs
	else
		echo -e "Cloning tlibs..."
		git clone https://code.ill.fr/scientific-software/takin/tlibs.git
	fi
fi


if [ -L takin2 ] || [ -d takin2 ]; then
	echo -e "A takin2 directory already exists. Skipping.";
else
	if [ -d archived/takin2 ] && [ $use_archived_copies -eq 1 ]; then
		echo -e "Using local archived copy of takin2.";
		ln -sf archived/takin2
	else
		echo -e "Cloning takin/mag-core repo..."
		git clone https://code.ill.fr/scientific-software/takin/mag-core.git
		mv -v mag-core takin2
	fi
fi


if [ -L tlibs2 ] || [ -d tlibs2 ]; then
	echo -e "A tlibs2 directory already exists. Skipping.";
else
	if [ -d archived/tlibs2 ] && [ $use_archived_copies -eq 1 ]; then
		echo -e "Using local archived copy of tlibs2.";
		ln -sf archived/tlibs2
	else
		echo -e "Cloning tlibs2 repo..."
		git clone https://code.ill.fr/scientific-software/takin/tlibs2.git
	fi
fi
