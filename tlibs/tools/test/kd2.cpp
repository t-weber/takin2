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

// test for k-d implementation: sort image into tree and recreate it again
// clang -march=native -O2 -o kd2 kd2.cpp -std=c++11 -lstdc++ -ljpeg

#include <iostream>
#include <vector>
#include <list>
#include <typeinfo>
#include <boost/gil/gil_all.hpp>
//#include <boost/gil/extension/io/png_io.hpp>
#include <boost/gil/extension/io/jpeg_dynamic_io.hpp>
#include "../math/kd.h"
#include "../log/debug.h"

namespace gil = boost::gil;

int main()
{
	gil::image<gil::rgb8_image_t, true> img;

	try
	{
		std::cout << "Loading image." << std::endl;
		gil::jpeg_read_image("/home/tweber/Pictures/rac.jpg", img);
	}
	catch(const std::ios_base::failure& ex)
	{
		std::cerr << "Error: Cannot load image" << std::endl;
		return 0;
	}

	auto view_gray = gil::color_converted_view<gil::gray8_pixel_t>(img._view);
	//std::cout << typeid(view_gray).name() << std::endl;

	gil::gray8s_pixel_t *pcPix = new gil::gray8s_pixel_t[view_gray.height() * view_gray.width()];
	gil::gray8s_view_t view = gil::interleaved_view(view_gray.width(), view_gray.height(), pcPix, view_gray.width());



	std::list<std::vector<int>> lst;


	for(int iY=0; iY<view_gray.height(); ++iY)
		for(int iX=0; iX<view_gray.width(); ++iX)
		{
			// *(view.row_begin(iY) + iX) = *(view_gray.row_begin(iY) + iX) - 128;
			lst.push_back({int(view_gray.width()-iX-1), iY, *(view_gray.row_begin(iY) + iX) - 128});
		}



	std::cout << "Generating k-d." << std::endl;
	tl::Kd<int> kd;
	kd.Load(lst, 2);


	std::cout << "Recreating image from k-d." << std::endl;
	for(int iY=0; iY<view_gray.height(); ++iY)
	{
		std::cout << "\rRow " << iY+1 << " / " << view_gray.height() << "." << std::flush;
		for(int iX=0; iX<view_gray.width(); ++iX)
		{
			const std::vector<int>& vecN = kd.GetNearestNode(std::vector<int>{iX, iY});
			*(view.row_begin(iY) + iX) = char(vecN[2]);
		}
	}
	std::cout << std::endl;


	std::cout << "Writing image." << std::endl;
	auto viewGray = gil::color_converted_view<gil::gray8_pixel_t>(view);
	std::cout << tl::get_typename<decltype(viewGray)>() << std::endl;

	gil::jpeg_write_view("/tmp/tst.jpg", viewGray, 85);
	delete[] pcPix;
}
