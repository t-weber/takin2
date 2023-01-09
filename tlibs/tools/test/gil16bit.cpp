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

// gcc -I . -o gil16bit gil16bit.cpp ../log/log.cpp -std=c++11 -lstdc++ -lpng -ltiff

#define NO_JPEG
#include "../gfx/gil.h"
#include "../log/debug.h"
#include <iostream>

namespace gil = tl::gil;

int main(int argc, char** argv)
{
	if(argc <= 1)
	{
		tl::log_err("Usage: ", argv[0], " <image>");
		return -1;
	}

	std::string strImg = argv[1];

	auto pimg = tl::load_image<gil::gray16_pixel_t, false>(strImg.c_str());
	if(!pimg)
		return -1;

	auto view = gil::view(*pimg);

	std::size_t iChan = view.num_channels();
	std::size_t iW = view.width();
	std::size_t iH = view.height();
	std::cout << iW << " x " << iH << ", " << iChan << " channels." << std::endl;

	//const auto iterRow = view.row_begin(908);
	//for(unsigned int iX=0; iX<iW; ++iX)
	//	std::cout << iterRow[iX] << " ";

	unsigned long lBkg = tl::get_roi_avg(view, 0, 500, 0, 300);
	std::cout << "Background counts: " << lBkg << std::endl;

	unsigned long lInt = tl::get_roi_sum(view, 310, 365, 880, 940, lBkg);
	std::cout << "Integrated counts: " << lInt << std::endl;

	return 0;
}
