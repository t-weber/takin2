#!/bin/bash
#
# @date 16-dec-2020
# @author Tobias Weber <tweber@ill.fr>
# @license GPLv2
#
# links to the system python libraries instead and removes the packaged ones
#

PRG="takin.app"

TOOL=install_name_tool
STRIP=strip

PY_VER=3.9


# files whose linkage is to be changed
declare -a filestochange=(
	"${PRG}/Contents/MacOS/takinmod_py"
)



# original symbols
declare -a changefrom=(
	"@executable_path/../Frameworks/Python.framework/Versions/${PY_VER}/Python"
)


# symbols to change into
declare -a changeto=(
	"/Library/Frameworks/Python.framework/Versions/Current/Python"
)


CNT=$(expr ${#changefrom[*]} - 1)


function fix_name()
{
	cfile=$1

	if [ ! -e ${cfile} ]; then
		echo -e "Warning: ${cfile} does not exist."
		return
	fi

	echo -e "Processing binary \"${cfile}\"..."
	chmod a+rx ${cfile}

	for idx in $(seq 0 ${CNT}); do
		cfrom=${changefrom[$idx]}
		cto=${changeto[$idx]}

		echo -e "\tChanging \"${cfrom}\"\n\t -> \"${cto}\"."
		chmod u+w ${cfile}
		${TOOL} -change ${cfrom} ${cto} ${cfile}
	done

	${STRIP} ${cfile}
	echo -e ""
}


# fix names for given files
for cfile in ${filestochange[@]}; do
	fix_name $cfile
done


# remove the packaged python frameworks and now non-needed libraries
rm -rfv "${PRG}/Contents/Frameworks/Python.framework"
rm -rfv "${PRG}/Contents/site-packages"
rm -fv "${PRG}/Contents/Libraries/libgfortran.5.dylib"
rm -fv "${PRG}/Contents/Libraries/libquadmath.0.dylib"
rm -fv "${PRG}/Contents/Libraries/libgomp.1.dylib"
rm -fv "${PRG}/Contents/Libraries/libgcc_s.1.dylib"
rm -fv "${PRG}/Contents/Libraries/libopenblas.0.dylib"
