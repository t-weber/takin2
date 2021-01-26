#!/bin/bash
#
# install all packages needed for building on focal (and bionic)
# @author Tobias Weber <tobias.weber@tum.de>
# @date 28-jul-20
# @license GPLv2
#

#if [[ $(id -u) -gt 0 ]]; then
#	echo -e "Please run this script as root."
#	exit -1
#fi


# -----------------------------------------------------------------------------
# install packages
# -----------------------------------------------------------------------------
if ! sudo apt-get install cmake clang build-essential \
	libboost-all-dev libclipper-dev \
	qt5-default qttools5-dev-tools libqt5svg5-dev qt5-assistant \
	libqwt-qt5-dev libpython3-dev \
	libfreetype6-dev libbz2-dev wget coreutils
then
	echo -e "Error: Could not install packages necessary for building."
	exit -1
fi
# -----------------------------------------------------------------------------
