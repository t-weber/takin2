/**
 * removes uneven rows in data files
 * @author Tobias Weber <tobias.weber@tum.de>
 * @license GPLv2
 * gcc -o clean_dat clean_dat.cpp -lstdc++ -std=c++11
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

#include <iostream>
#include <fstream>
#include <string>
#include <tuple>
#include "tlibs/string/string.h"


std::tuple<bool, std::size_t> cleanfile(const char* pcFile)
{
	std::ifstream ifstr(pcFile);
	if(!ifstr)
	{
		std::cerr << "Could not open " << pcFile << " for reading." << std::endl;
		return std::make_tuple(false,0);
	}

	const int iMinCols = 3;
	const bool bExactColNumMatch = 0;
	const bool bShrinkToMatchCols = 0;

	int iCols = -1;
	std::vector<std::vector<std::string>> vecRows;

	std::size_t iLinesRemoved = 0;
	std::size_t iLine = 1;
	while(1)
	{
		std::string strLine;
		std::getline(ifstr, strLine);
		if(ifstr.eof())
			break;

		std::vector<std::string> vecToks;
		tl::get_tokens<std::string>(strLine, std::string(" \t"), vecToks);

		if(iCols < 0)
		{
			iCols = int(vecToks.size());
			//std::cout << "First line has " << iCols << " columns." << std::endl;
		}

		if(bShrinkToMatchCols && vecToks.size() > iCols)
			vecToks.resize(iCols);

		bool bInsert = 0;
		if(vecToks.size() > iMinCols)
		{
			if(bExactColNumMatch && vecToks.size() == iCols)
				bInsert = 1;
			if(!bExactColNumMatch)
				bInsert = 1;

			if(bInsert)
				vecRows.emplace_back(std::move(vecToks));
		}

		if(!bInsert)
		{
			//std::cout << "Removing line " << iLine << "." << std::endl;
			++iLinesRemoved;
		}

		++iLine;
	}
	ifstr.close();


	std::ofstream ofstr(pcFile);
	if(!ofstr)
	{
		std::cerr << "Could not open " << pcFile << " for writing." << std::endl;
		return std::make_tuple(false,0);
	}

	for(const std::vector<std::string>& vecRow : vecRows)
	{
		for(const std::string& str : vecRow)
			ofstr << str << " ";
		ofstr << "\n";
	}

	return std::make_tuple(true, iLinesRemoved);
}


int main(int argc, char** argv)
{
	if(argc < 2)
	{
		std::cerr << "Usage: " << argv[0] << " <file1.dat> <file2.dat> ..." << std::endl;
		return false;
	}

	for(int iArg=1; iArg<argc; ++iArg)
	{
		std::cout << argv[iArg] << " ... ";

		bool bOk = 0;
		std::size_t iLinesRemoved = 0;
		std::tie(bOk, iLinesRemoved) = cleanfile(argv[iArg]);

		if(bOk)
			std::cout << " OK, " << iLinesRemoved << " lines removed."  << std::endl;
		else
			std::cout << " ERROR." << std::endl;
	}

	return 0;
}
