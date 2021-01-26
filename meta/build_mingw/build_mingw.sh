#!/bin/bash
#
# takin mingw build script
# @author Tobias Weber <tweber@ill.fr>
# @date sep-2020
# @license GPLv2
#

# individual building steps
setup_buildenv=0
setup_externals=1
setup_externals2=1
build_externals=1
build_takin=0
build_takin2=0
build_package=1


NUM_CORES=$(nproc)


# get root dir of takin repos
TAKIN_ROOT=$(dirname $0)/../..
cd "${TAKIN_ROOT}"
TAKIN_ROOT=$(pwd)
echo -e "Takin root dir: ${TAKIN_ROOT}"


if [ $setup_buildenv -ne 0 ]; then
	echo -e "\n================================================================================"
	echo -e "Setting up build environment..."
	echo -e "================================================================================\n"

	pushd "${TAKIN_ROOT}/core"
	# TODO: set up mingw packages
	#	./setup_lin/buildenv_focal.sh
	popd
fi


if [ $setup_externals -ne 0 ]; then
	echo -e "\n================================================================================"
	echo -e "Getting external dependencies (1/2)..."
	echo -e "================================================================================\n"

	pushd "${TAKIN_ROOT}/core"
		rm -rf tmp
		./setup/setup_externals.sh
	popd
fi


if [ $setup_externals2 -ne 0 ]; then
	echo -e "\n================================================================================"
	echo -e "Getting external dependencies (2/2)..."
	echo -e "================================================================================\n"

	pushd "${TAKIN_ROOT}/mag-core"
		rm -rf ext
		./setup/setup_externals.sh
	popd
fi


if [ $build_externals -ne 0 ]; then
	echo -e "\n================================================================================"
	echo -e "Building external libraries (Minuit)..."
	echo -e "================================================================================\n"

	pushd "${TAKIN_ROOT}/tmp"
		"${TAKIN_ROOT}"/meta/externals/build_minuit.sh --mingw
	popd
fi


if [ $build_takin -ne 0 ]; then
	echo -e "\n================================================================================"
	echo -e "Building main Takin binary..."
	echo -e "================================================================================\n"

	pushd "${TAKIN_ROOT}/core"
		./setup/clean.sh

		mkdir -p build
		cd build
		mingw64-cmake -DDEBUG=False ..
		mingw64-make #-j${NUM_CORES}
	popd
fi


if [ $build_takin2 -ne 0 ]; then
	echo -e "\n================================================================================"
	echo -e "Building Takin 2 tools..."
	echo -e "================================================================================\n"

	pushd "${TAKIN_ROOT}/mag-core"
		rm -rf build
		mkdir -p build
		cd build

		mingw64-cmake -DCMAKE_BUILD_TYPE=Release -DONLY_BUILD_FINISHED=True ..
		mingw64-make #-j${NUM_CORES}

		# copy tools to Takin main dir
		cp -v tools/cif2xml/takin_cif2xml.exe "${TAKIN_ROOT}"/core/bin/
		cp -v tools/cif2xml/takin_findsg.exe "${TAKIN_ROOT}"/core/bin/
		cp -v tools/pol/takin_pol.exe "${TAKIN_ROOT}"/core/bin/
	popd
fi


if [ $build_package -ne 0 ]; then
	echo -e "\n================================================================================"
	echo -e "Building Takin package..."
	echo -e "================================================================================\n"

	pushd "${TAKIN_ROOT}"
		rm -rf tmp-mingw
		cd core
		./setup_mingw/cp_mingw_takin.sh "${TAKIN_ROOT}/tmp-mingw/takin"

		cd ../tmp-mingw
		zip -9 -r takin.zip takin
	popd


	if [ -e  "${TAKIN_ROOT}/tmp-mingw/takin.zip" ]; then
		echo -e "\n================================================================================"
		echo -e "The built Takin package can be found here:\n\t${TAKIN_ROOT}/tmp/takin.zip"
		echo -e "================================================================================\n"
	else
		echo -e "\n================================================================================"
		echo -e "Error: Takin package could not be built!"
		echo -e "================================================================================\n"
	fi
fi
