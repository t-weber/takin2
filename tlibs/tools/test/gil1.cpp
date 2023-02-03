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

// gcc -I . -o gil1 gil1.cpp ../log/log.cpp -std=c++11 -lstdc++ -ljpeg -lpng -ltiff

#include "../gfx/gil.h"
#include "../log/debug.h"
#include <iostream>

namespace gil = tl::gil;

int main()
{
	auto pimg = tl::load_image<gil::rgb8_image_t>("/home/tweber/Pictures/rac.jpg");
	if(!pimg)
		return -1;

	std::size_t iChan = gil::view(*pimg).num_channels();
	std::size_t iW = gil::view(*pimg).width();
	std::size_t iH = gil::view(*pimg).height();
	std::cout << iW << " x " << iH << ", " << iChan << " channels." << std::endl;

	tl::for_each_in_view(gil::view(*pimg),
		[](gil::rgb8_pixel_t& pix)
		{
			std::swap(pix[0], pix[2]);
		});

	if(!tl::save_image("/tmp/tst.png", pimg/*.get()*/))
		return -1;

	return 0;
}
