/**
 * tlibs test file
 * @author Tobias Weber <tweber@ill.fr>
 * @date 18/oct/2022
 * @license GPLv2 or GPLv3
 *
 * ----------------------------------------------------------------------------
 * tlibs -- a physical-mathematical C++ template library
 * Copyright (C) 2017-2022  Tobias WEBER (Institut Laue-Langevin (ILL),
 *                          Grenoble, France).
 * Copyright (C) 2015-2017  Tobias WEBER (Technische Universitaet Muenchen
 *                          (TUM), Garching, Germany).
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 2 of the License, or
 * (at your option) version 3.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 * ----------------------------------------------------------------------------
 */

// g++ -std=c++11 -I../.. -I /usr/include/hdf5/serial -o h5 h5.cpp -L /usr/lib/x86_64-linux-gnu/hdf5/serial/ -lhdf5_cpp -lhdf5

#include <iostream>
#include <string>
#include "file/h5.h"

int main()
{
	H5::H5File h5file("/users/tw/Downloads/mail_tmp/thales.nxs", H5F_ACC_RDONLY);

	std::vector<std::string> entries;
	tl::get_h5_entries(h5file, "/", entries);
	for(const std::string& entry : entries)
	{
		std::string nx_class = tl::get_h5_attr<std::string>(h5file, "/" + entry, "NX_class");
		std::cout << entry << " -- " << nx_class << std::endl;
	}
	std::cout << std::endl;

	std::cout << std::endl;

	int steps = 0;
	bool ok = tl::get_h5_scalar(h5file, "entry0/data_scan/total_steps", steps);
	std::cout << std::boolalpha << ok << std::endl;
	std::cout << "steps: " << steps << std::endl;

	std::cout << std::endl;

	std::string ident;
	ok = tl::get_h5_string(h5file, "entry0/experiment_identifier", ident);
	std::cout << std::boolalpha << ok << std::endl;
	std::cout << "ident: " << ident << std::endl;

	std::cout << std::endl;

	std::vector<int> scanned;
	ok = tl::get_h5_vector(h5file, "entry0/data_scan/scanned_variables/variables_names/scanned", scanned);
	std::cout << std::boolalpha << ok << std::endl;
	for(int i : scanned)
	std::cout << i << " ";
	std::cout << std::endl;

	std::cout << std::endl;

	std::vector<std::string> property;
	ok = tl::get_h5_string_vector(h5file, "entry0/data_scan/scanned_variables/variables_names/property", property);
	std::cout << std::boolalpha << ok << std::endl;
	for(const std::string& str : property)
	std::cout << str << " ";
	std::cout << std::endl;

	std::cout << std::endl;

	std::vector<std::vector<double>> data;
	ok = tl::get_h5_matrix(h5file, "entry0/data_scan/scanned_variables/data", data);
	std::cout << std::boolalpha << ok << std::endl;
	for(const auto& row : data)
	{
		for(double val : row)
			std::cout << val << " ";
		std::cout << std::endl;
	}

	h5file.close();
	return 0;
}
