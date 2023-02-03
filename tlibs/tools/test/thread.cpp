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

// clang -o thread thread.cpp -std=c++11 -lstdc++ -lpthread

#include "../helper/thread.h"
#include <iostream>

int tst(int i)
{
	return i*i;
}

int main()
{
	tl::ThreadPool<int()> tp;
	for(int i=0; i<10; ++i)
		tp.AddTask(std::bind(tst, i));

	for(int i=0; i<10; ++i)
		tp.AddTask([i]()->int{ return i; });

	tp.StartTasks();

	auto& lstFut = tp.GetFutures();
	for(auto& fut : lstFut)
		std::cout << fut.get() << std::endl;
	return 0;
}
