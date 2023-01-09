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

// gcc -o proctst proctst.cpp -lstdc++ -lboost_iostreams

#include "../helper/proc.h"
#include "../string/string.h"
#include "../time/chrono.h"
#include <iostream>
#include <fstream>
#include <thread>

using t_real = double;

int main(int argc, char **argv)
{
	std::ofstream ofstr("sensors.dat");
	std::string strStartTime = tl::epoch_to_str(tl::epoch<t_real>());
	ofstr << "#\n# start time: " << strStartTime << "\n#\n";

	std::size_t iRun = 0;
	while(1)
	{
		std::cout << "\r" << (iRun+1);
		std::cout.flush();

		ofstr << iRun << " \t";
		tl::PipeProc<> proc("sensors", 0);
		if(proc.IsReady())
		{
			std::istream& istr = proc.GetIstr();
			while(!istr.eof())
			{
				std::string strLine;
				std::getline(istr, strLine);

				if(strLine.find("crit") != std::string::npos)
				{
					std::string str = tl::str_between(strLine, std::string(":"), std::string("("), 1, 0);

					t_real dVar = tl::str_to_var<t_real>(str);
					ofstr << dVar << " \t";
				}
			}
		}

		ofstr << std::endl;
		std::this_thread::sleep_for(std::chrono::milliseconds(1000));

		++iRun;
	}

	return 0;
}
