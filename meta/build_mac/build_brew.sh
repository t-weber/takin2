#!/bin/bash
#
# takin brew build script
# @author Tobias Weber <tweber@ill.fr>
# @date sep-2020
# @license GPLv2
#


# individual building steps
setup_buildenv=1
setup_externals=1
setup_externals2=1
build_takin=1
build_takin2=1
build_package=1

use_syspy=1


export MACOSX_DEPLOYMENT_TARGET=10.10

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
		./setup_mac/buildenv_brew.sh
	popd
fi


if [ $setup_externals -ne 0 ]; then
	echo -e "\n================================================================================"
	echo -e "Getting external dependencies (1/2)..."
	echo -e "================================================================================\n"

	pushd "${TAKIN_ROOT}/core"
		rm -rf tmp
		./setup/setup_externals.sh
		./setup/get_3rdparty_licenses.sh
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


if [ $build_takin -ne 0 ]; then
	echo -e "\n================================================================================"
	echo -e "Building main Takin binary..."
	echo -e "================================================================================\n"

	pushd "${TAKIN_ROOT}/core"
		./setup/clean.sh

		mkdir -p build
		cd build
		cmake -DDEBUG=False ..
		make -j${NUM_CORES}
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
		cmake -DCMAKE_BUILD_TYPE=Release -DONLY_BUILD_FINISHED=True ..
		make -j${NUM_CORES}

		# copy tools to Takin main dir
		cp -v tools/cif2xml/takin_cif2xml "${TAKIN_ROOT}"/core/bin/
		cp -v tools/cif2xml/takin_findsg "${TAKIN_ROOT}"/core/bin/
		cp -v tools/pol/takin_pol "${TAKIN_ROOT}"/core/bin/
		cp -v tools/structfact/takin_structfact "${TAKIN_ROOT}"/core/bin/
		cp -v tools/magstructfact/takin_magstructfact "${TAKIN_ROOT}"/core/bin/
		cp -v tools/scanbrowser/takin_scanbrowser "${TAKIN_ROOT}"/core/bin/
		cp -v tools/magsgbrowser/takin_magsgbrowser "${TAKIN_ROOT}"/core/bin/
		cp -v tools/moldyn/takin_moldyn "${TAKIN_ROOT}"/core/bin/
	popd
fi


if [ $build_package -ne 0 ]; then
	pushd "${TAKIN_ROOT}"
		rm -rf tmp
		mkdir -p tmp

		cd core

		echo -e "\n================================================================================"
		echo -e "Creating icons..."
		echo -e "================================================================================\n"
		./setup_mac/mk_icon.sh

		echo -e "\n================================================================================"
		echo -e "Creating app directory..."
		echo -e "================================================================================\n"
		./setup_mac/cp_app.sh

		echo -e "\n================================================================================"
		echo -e "Fixing dynamic binding for local libraries..."
		echo -e "================================================================================\n"
		./setup_mac/fix_names.sh

		if [ $use_syspy -ne 0 ]; then
			echo -e "\n================================================================================"
			echo -e "Using system python frameworks instead..."
			echo -e "================================================================================\n"
			./setup_mac/use_syspy.sh
		fi

		echo -e "\n================================================================================"
		echo -e "Building Takin package..."
		echo -e "================================================================================\n"
		./setup_mac/cp_dmg.sh

		cp -v takin.dmg ../tmp/takin.dmg
	popd


	if [ -e  "${TAKIN_ROOT}/tmp/takin.dmg" ]; then
		echo -e "\n================================================================================"
		echo -e "The built Takin package can be found here:\n\t${TAKIN_ROOT}/tmp/takin.dmg"
		echo -e "================================================================================\n"
	else
		echo -e "\n================================================================================"
		echo -e "Error: Takin package could not be built!"
		echo -e "================================================================================\n"
	fi
fi

