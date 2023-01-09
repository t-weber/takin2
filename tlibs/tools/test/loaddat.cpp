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

// clang++ -DNO_IOSTR -o loaddat test/loaddat.cpp ../log/log.cpp -std=c++11 -lboost_iostreams -lm

#include <iostream>
#include "../file/loaddat.h"

int main()
{
	//tl::comp_stream_to_stream<char>(std::cin, std::cout);
	//tl::comp_stream_to_stream<wchar_t>(std::wcin, std::wcout);

	//std::wstring str(L"123");
	//std::wcout << tl::str_to_var<int, std::wstring>(str) << std::endl;

	//tl::DatFile<float, wchar_t> dat;
	tl::DatFile<float, char> dat;

	dat.SetSeparatorChars({':'});
	//dat.Load("/home/tw/Measurements/mira-mnsi-14/data/9925_00009990.dat");
	//dat.Load(L"/home/tw/Measurements/mira-mnsi-14/data/9925_00009990.dat");
	dat.Load("/home/tweber/Messdaten/mira-mnsi-15/data/10365_00015382.dat");

	if(dat.IsOk())
	{
		std::cout << "Number of columns: " << dat.GetColumnCount() << std::endl;
		std::cout << "Number of rows: " << dat.GetColumn(0).size() << std::endl;

		const auto& hdr = dat.GetHeader();
		for(const auto& pair : hdr)
		{
			std::cout << pair.first << " = " << pair.second << std::endl;
			//std::wcout << pair.first << L" = " << pair.second << std::endl;
		}

		dat.Save("/tmp/tst.dat");
	}

	return 0;
}
