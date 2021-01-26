/**
 * metropolis test
 * @author Tobias Weber <tobias.weber@tum.de>
 * @date 8-oct-16
 * @license GPLv2 or GPLv3
 */
// g++ -march=native -O2 -std=c++11 -DNO_JPEG -DNO_TIFF -o metrop metrop.cpp ../../log/log.cpp ../../math/rand.cpp -lpng

#include "../../phys/mag.h"
#include "../../math/rand.h"
#include "../../phys/units.h"
#include "../../gfx/gil.h"
#include "../../helper/array.h"
#include <iostream>

using t_real = double;
static const tl::t_energy_si<t_real> meV = tl::get_one_meV<t_real>();
static const tl::t_temperature_si<t_real> kelvin = tl::get_one_kelvin<t_real>();
static const t_real k = t_real(tl::get_kB<t_real>() / meV * kelvin);

int main()
{
	tl::init_rand();

	long iW = 256;
	long iH = 256;

	t_real J = 1.;		// meV
	t_real T = 25.;		// K

	std::size_t iIter = 1e7;

	std::cout << "Width: "; std::cin >> iW;
	std::cout << "Height: "; std::cin >> iH;
	std::cout << "Iterations: "; std::cin >> iIter;
	std::cout << "J (meV): "; std::cin >> J;
	std::cout << "T (K): "; std::cin >> T;

	boost::multi_array<bool, 2> arr
		= tl::rand_array<bool, 2, boost::array, boost::multi_array>({iW, iH});

	t_real Etot = 0.;
	tl::log_info("Running Metropolis algo...");
	tl::metrop<t_real, 2>({iW,iH}, iIter, J, k, T, arr, &Etot);
	tl::log_info("Finished Metropolis algo. E/N = ", Etot, " meV.");


	tl::log_info("Writing image.");

	using t_view = tl::gil::gray8_view_t;
	using t_pix = typename t_view::value_type;
	t_view view;
	std::vector<t_pix> vecPix;
	tl::create_imgview<t_view>(iW, iH, vecPix, view);

	for(long i=0; i<iH; ++i)
		for(long j=0; j<iW; ++j)
			*(view.row_begin(i)+j) = arr[i][j]*255;

	tl::save_view("metrop.png", &view);
	return 0;
}
