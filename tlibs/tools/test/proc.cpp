/**
 * tlibs test file
 * @author Tobias Weber <tobias.weber@tum.de>
 * @license GPLv2 or GPLv3
 *
 * ----------------------------------------------------------------------------
 * tlibs -- a physical-mathematical C++ template library
 * Copyright (C) 2017-2021  Tobias WEBER (Institut Laue-Langevin (ILL),
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

// gcc -o proc proc.cpp -lstdc++ -lboost_iostreams

#include "../helper/proc.h"
#include <iostream>

int main(int argc, char **argv)
{
	if(argc < 2) return -1;

	tl::PipeProc<> proc(argv[1], 0);
	if(proc.IsReady())
	{
		std::istream& istr = proc.GetIstr();
		while(!istr.eof())
		{
			std::string strLine;
			std::getline(istr, strLine);
			std::cout << strLine << std::endl;
		}
	}

	return 0;
}
