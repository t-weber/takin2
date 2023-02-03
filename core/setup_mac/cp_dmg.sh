#!/bin/bash
#
# @date jan-19
# @author Tobias Weber <tobias.weber@tum.de>
# @license GPLv2
#
# copy application into disk image
#
# ----------------------------------------------------------------------------
# Takin (inelastic neutron scattering software package)
# Copyright (C) 2017-2023  Tobias WEBER (Institut Laue-Langevin (ILL),
#                          Grenoble, France).
# Copyright (C) 2013-2017  Tobias WEBER (Technische Universitaet Muenchen
#                          (TUM), Garching, Germany).
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; version 2 of the License.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along
# with this program; if not, write to the Free Software Foundation, Inc.,
# 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
# ----------------------------------------------------------------------------
#

echo -e "Removing any old file ..."
rm -f takin.dmg

echo -e "Cleaning temporary files ..."
find takin.app -name ".DS_Store" -exec rm -fv {} \;
find takin.app -name ".dir" -exec rm -fv {} \;
find takin.app -type d -name "__pycache__" -exec rm -rfv {} \;

echo -e "\nCreating a writable image ..."
if ! hdiutil create takin.dmg -srcfolder takin.app -fs UDF -format "UDRW" -volname "takin"; then
	echo -e "Cannot create image file."
	exit -1
fi

echo -e "\n--------------------------------------------------------------------------------"

echo -e "Mounting the image ..."
if ! hdiutil attach takin.dmg -readwrite; then
	echo -e "Cannot mount image file."
	exit -1
fi

echo -e "\n--------------------------------------------------------------------------------"

echo -e "Creating additional files/links in the image ..."
ln -sf /Applications /Volumes/takin/Install_by_dragging_Takin_here

echo -e "\n--------------------------------------------------------------------------------"

echo -e "Unmounting the image ..."
if ! hdiutil detach /Volumes/takin; then
	echo -e "Cannot detach image file."
	exit -1
fi

echo -e "\n--------------------------------------------------------------------------------"

echo -e "Converting the image to read-only and compressing ..."
if ! hdiutil convert takin.dmg -o takin_final.dmg -format "UDBZ"; then
	echo -e "Cannot convert image file to UDBZ."
	exit -1
fi

echo -e "\n--------------------------------------------------------------------------------"

echo -e "Copying file ..."
mv -v takin_final.dmg takin.dmg

echo -e "\nSuccessfully created takin.dmg."
