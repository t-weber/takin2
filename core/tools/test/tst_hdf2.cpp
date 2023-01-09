/**
 * hdf test
 * @author Tobias Weber <tweber@ill.fr>
 * @date 12-nov-21
 * @license GPLv2
 *
 * Reference: https://support.hdfgroup.org/HDF5/doc1.6/UG/11_Datatypes.html
 *
 * g++ -o tst_hdf2 tst_hdf2.cpp -lhdf5_serial
 *
 * ----------------------------------------------------------------------------
 * Takin (inelastic neutron scattering software package)
 * Copyright (C) 2017-2021  Tobias WEBER (Institut Laue-Langevin (ILL),
 *                          Grenoble, France).
 * Copyright (C) 2013-2017  Tobias WEBER (Technische Universitaet Muenchen
 *                          (TUM), Garching, Germany).
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; version 2 of the License.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along
 * with this program; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 * ----------------------------------------------------------------------------
 */

#include <hdf5/serial/hdf5.h>
#include <iostream>
#include <vector>


int main(int argc, char **argv)
{
	if(argc <= 1)
	{
		std::cerr << "Please specify a nexus/hdf5 file." << std::endl;
		return -1;
	}

	hid_t file = ::H5Fopen(argv[1], H5F_ACC_RDONLY, H5P_DEFAULT);
	if(file < 0)
	{
		std::cerr << "Cannot open file \"" << argv[1] << "\"." << std::endl;
		return -1;
	}

	hid_t data = ::H5Dopen(file, "/entry0/data_scan/scanned_variables/variables_names/name", H5P_DEFAULT);
	if(data < 0)
	{
		std::cerr << "Cannot open data." << std::endl;
		return -1;
	}
	hid_t set = ::H5Dget_space(data);
	if(set < 0)
	{
		std::cerr << "Cannot create data set." << std::endl;
		return -1;
	}

	int iRank = ::H5Sget_simple_extent_ndims(set);
	std::vector<hsize_t> vecDim(iRank);
	::H5Sget_simple_extent_dims(set, vecDim.data(), 0);
	::H5Sclose(set);
	std::cout << "Rank: " << iRank << ", dims: ";
	for(hsize_t iDim : vecDim)
		std::cout << iDim << ", ";
	std::cout << std::endl;

	hid_t type_id = ::H5Tcopy(H5T_C_S1);
	::H5Tset_size(type_id, H5T_VARIABLE);

	hsize_t num_strings = vecDim[0];
	//hvl_t *str = new hvl_t[num_strings];
	const char **str = new const char*[num_strings];
	if(::H5Dread(data, type_id, H5S_ALL, H5S_ALL, H5P_DEFAULT, str) < 0)
	{
		std::cerr << "Cannot read data." << std::endl;
		return -1;
	}

	for(hsize_t i=0; i<num_strings; ++i)
	{
		std::cout << "String " << i << ": " << str[i] << std::endl;
	}

	delete[] str;

	::H5Dclose(data);
	::H5Fclose(file);

	::H5close();
	return 0;
}
