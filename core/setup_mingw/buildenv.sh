#!/bin/bash
#
# install all packages needed for building on mingw under fedora
# @author Tobias Weber <tweber@ill.fr>
# @date jan-2021
# @license GPLv2
#

if [[ $(id -u) -gt 0 ]]; then
	echo -e "Please run this script as root."
	exit -1
fi


# -----------------------------------------------------------------------------
# install packages
# -----------------------------------------------------------------------------
if ! dnf install mingw64-gcc mingw64-gcc-c++ mingw64-boost \
	mingw64-qt5-qtbase mingw64-qt5-qtbase-devel mingw64-qt5-qttools \
	mingw64-qt5-qttools-tools mingw64-qwt-qt5 \
	mingw64-python3 mingw64-python3-numpy \
	mingw64-freetype
then
	echo -e "Error: Could not install packages necessary for building."
	exit -1
fi
# -----------------------------------------------------------------------------
