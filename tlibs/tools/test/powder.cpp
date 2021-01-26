/**
 * tlibs test file
 * @author Tobias Weber <tobias.weber@tum.de>
 * @license GPLv2 or GPLv3
 *
 * gcc -o powder powder.cpp ../log/log.cpp -lstdc++ -lm -lMinuit2 -lgomp
 */

#include "../phys/powder.h"
#include <iostream>

int main()
{
	tl::Powder<int> powder;

	powder.AddPeak(1,1,0);
	powder.AddPeak(1,-1,0);
	powder.AddPeak(-1,-1,0);

	powder.AddPeak(1,0,0);
	powder.AddPeak(0,1,0);

	powder.AddPeak(2,1,0);
	powder.AddPeak(-4,0,0);

	std::cout << powder << std::endl;
	std::cout << std::endl;


	std::vector<double> vecGs = {1.8, 1.3, 0.9};
	std::vector<double> vecTTs = {78./180.*M_PI, 55./180.*M_PI, 40./180.*M_PI};
	std::vector<double> vecTTErrs = {5., 5., 5.};
	std::vector<double> vecRes = {1.4, 0., 0.};
	std::vector<double> vecResErr = {0.1, M_PI, 0.};
	tl::powder_align(3.4, vecGs, vecTTs, vecTTErrs, vecRes, vecResErr);

	std::cout << "mono twotheta: " << vecRes[2]/M_PI*180. << std::endl;
	return 0;
}
