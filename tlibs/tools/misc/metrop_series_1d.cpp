/**
 * 1d metropolis test (without phase transition)
 * @author Tobias Weber <tobias.weber@tum.de>
 * @date 8-oct-16
 * @license GPLv2 or GPLv3
 */
// g++ -march=native -O2 -std=c++11 -o metrop_series_1d metrop_series_1d.cpp ../../log/log.cpp ../../math/rand.cpp

#include "../../math/rand.h"
#include "../../phys/mag.h"
#include "../../phys/units.h"
#include "../../helper/array.h"
#include <iostream>
#include <fstream>

using t_real = double;
static const tl::t_energy_si<t_real> meV = tl::get_one_meV<t_real>();
static const tl::t_temperature_si<t_real> kelvin = tl::get_one_kelvin<t_real>();
static const t_real k = t_real(tl::get_kB<t_real>() / meV * kelvin);


int main()
{
	tl::init_rand();

	std::size_t iIter = 1e7;
	long iW = 256;
	t_real J = 2.5;		// meV

	boost::multi_array<bool, 1> arr
		= tl::rand_array<bool, 1, boost::array, boost::multi_array>({ iW });
	std::ofstream ofstrE("E1d.dat");

	ofstrE << std::setw(16) << std::left << "#T" << " " 
		<< std::setw(16) << std::left << "E_t" << " " 
		<< std::setw(16) << std::left << "<S>" << std::endl;

	for(t_real T=200.; T>=0.; T-=10.)
	{
		t_real Etot = 0., S = 0.;
		tl::log_info("T = ", T, " K");
		tl::metrop<t_real, 1>({iW}, iIter, J, k, T, arr, &Etot, &S);
		ofstrE << std::setw(16) << std::left << T << " " 
			<< std::setw(16) << std::left << Etot << " " 
			<< std::setw(16) << std::left << S << std::endl;
	}

	return 0;
}
