#!/bin/bash
#
# gets license files for 3rd party libraries
# @author Tobias Weber <tweber@ill.fr>
# @date jan-2021
# @license GPLv2
#

LICDIR=3rdparty_licenses

echo -e "Downloading license texts...\n"

# boost
if ! wget http://www.boost.org/LICENSE_1_0.txt -O ${LICDIR}/boost_license.txt; then
	echo -e "Error: Cannot download Boost license.";
fi

# python
if ! wget https://raw.githubusercontent.com/python/cpython/master/Doc/license.rst -O ${LICDIR}/python_license.txt; then
	echo -e "Error: Cannot download Python license.";
fi

# numpy
if ! wget https://raw.githubusercontent.com/numpy/numpy/master/LICENSE.txt -O ${LICDIR}/numpy_license.txt; then
	echo -e "Error: Cannot download Numpy license.";
fi

# scipy
if ! wget https://raw.githubusercontent.com/scipy/scipy/master/LICENSE.txt -O ${LICDIR}/scipy_license.txt; then
	echo -e "Error: Cannot download Scipy license.";
fi

# ROOT
if ! wget https://raw.githubusercontent.com/root-project/root/master/LICENSE -O ${LICDIR}/ROOT_license.txt; then
	echo -e "Error: Cannot download ROOT license.";
fi

# qwt
if ! wget https://raw.githubusercontent.com/qwtplot/qwt/trunk/COPYING -O ${LICDIR}/qwt_license.txt; then
	echo -e "Error: Cannot download Qwt license.";
fi

# freetype
if ! wget https://git.savannah.gnu.org/cgit/freetype/freetype2.git/plain/docs/FTL.TXT -O ${LICDIR}/freetype_license.txt; then
	echo -e "Error: Cannot download Freetype license.";
fi

# lapack(e)
if ! wget http://www.netlib.org/lapack/LICENSE.txt -O ${LICDIR}/lapack_license.txt; then
	echo -e "Error: Cannot download Lapack(e) license.";
fi

# qhull
if ! wget https://raw.githubusercontent.com/qhull/qhull/master/COPYING.txt -O ${LICDIR}/qhull_license.txt; then
	echo -e "Error: Cannot download Qhull license.";
fi

# neutronpy
if ! wget https://raw.githubusercontent.com/neutronpy/neutronpy/master/LICENSE.txt -O ${LICDIR}/neutronpy_license.txt; then
	echo -e "Error: Cannot download NeutronPy license.";
fi

# bodr
if ! wget https://raw.githubusercontent.com/egonw/bodr/master/bodr/COPYING -O ${LICDIR}/blue_obelisk_data_repository_license.txt; then
	echo -e "Error: Cannot download BODR license.";
fi

# dejavu
if ! wget https://dejavu-fonts.github.io/License.html -O ${LICDIR}/dejavu_license.htm; then
	echo -e "Error: Cannot download DejaVu license.";
fi

# libjpeg
if ! wget https://raw.githubusercontent.com/freedesktop/libjpeg/master/README -O ${LICDIR}/libjpeg_license.txt; then
	echo -e "Error: Cannot download libjpg license.";
fi

# libtiff
if ! wget http://www.libtiff.org/misc.html -O ${LICDIR}/libtiff_license.htm; then
	echo -e "Error: Cannot download libtiff license.";
fi

# libpng
if ! wget http://www.libpng.org/pub/png/src/libpng-LICENSE.txt -O ${LICDIR}/libpng_license.txt; then
	echo -e "Error: Cannot download libpng license.";
fi

# gnuplot
if ! wget https://raw.githubusercontent.com/gnuplot/gnuplot/master/Copyright -O ${LICDIR}/gnuplot_license.txt; then
	echo -e "Error: Cannot download Gnuplot license.";
fi

# gemmi
if ! wget https://raw.githubusercontent.com/project-gemmi/gemmi/master/LICENSE.txt -O ${LICDIR}/gemmi_license.txt; then
	echo -e "Error: Cannot download Gemmi license.";
fi

# qcustomplot
if ! wget https://gitlab.com/DerManu/QCustomPlot/-/raw/master/GPL.txt -O ${LICDIR}/qcustomplot_license.txt; then
	echo -e "Error: Cannot download QCustomPlot license.";
fi
