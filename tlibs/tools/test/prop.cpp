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

// clang -o prop prop.cpp ../log/log.cpp -std=c++11 -lstdc++ -lm -lboost_iostreams

#include "../file/prop.h"

int main()
{
	tl::Prop<std::string> prop;
	//prop.Load("tst.ini");
	prop.Add<std::string>("Sec1/Key1", "123");
	prop.Add<std::string>("Sec1/Key2", "Test");

	std::map<std::string, std::string> map;
	map["Sec2/Key1"] = "abcde";
	map["Sec2/Key2"] = "xyz";
	prop.Add(map);

	prop.Save("tst.xml");

	std::cout << prop.Query<std::string>("Sec1/Key1") << std::endl;
	return 0;
}
