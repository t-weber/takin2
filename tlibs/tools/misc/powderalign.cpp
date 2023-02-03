/**
 * tlibs test file
 * @author Tobias Weber <tobias.weber@tum.de>
 * @license GPLv2 or GPLv3
 *
 * gcc -o powderalign powderalign.cpp ../../log/log.cpp -lstdc++ -lm -lMinuit2 -lgomp
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

#include "../../phys/powder.h"
#include "../../string/string.h"
#include <iostream>
#include <fstream>

using t_real = double;


int main()
{
	t_real dMono = 3.5, dKi = 1.4;
	std::string strd, strki;
	std::string strGs, strTTs, strTTErrs;

	std::cout << "Monochromator d [A] = ";
	std::getline(std::cin, strd);
	dMono = tl::str_to_var<t_real>(strd);

	std::cout << "Nominal ki [1/A] = ";
	std::getline(std::cin, strki);
	dKi = tl::str_to_var<t_real>(strki);

	std::cout << "Comma-separated theoretical powder line Gs [1/A] = ";
	std::getline(std::cin, strGs);

	std::cout << "Comma-separated measured powder line two-thetas [deg] = ";
	std::getline(std::cin, strTTs);

	std::cout << "Comma-separated measured powder line two-theta errors [deg] = ";
	std::getline(std::cin, strTTErrs);


	std::vector<t_real> vecGs, vecTTs, vecTTErrs;
	tl::get_tokens<t_real, std::string>(strGs, ",;", vecGs);
	tl::get_tokens<t_real, std::string>(strTTs, ",;", vecTTs);
	tl::get_tokens<t_real, std::string>(strTTErrs, ",;", vecTTErrs);

	for(t_real& dTT : vecTTs) dTT *= M_PI/180.;
	for(t_real& dTT : vecTTErrs) dTT *= M_PI/180.;

	if(vecGs.size() != vecTTs.size() || vecTTs.size() != vecTTErrs.size())
	{
		std::cerr << "Input array size mismatch!" << std::endl;
		return -1;
	}

	std::cout << vecGs.size() << " values entered." << std::endl;


	std::vector<t_real> vecRes = { dKi, 0., 0. },
		vecResErr = { dKi*0.25, 0.5, vecRes[2]*0.25 };

	tl::powder_align(dMono, vecGs, vecTTs, vecTTErrs, vecRes, vecResErr);


	std::cout << "Monochromator ki = " << vecRes[0] << "/A" << " (nominal: " << dKi << "/A)\n";
	std::cout << "Monochromator two-theta = " << vecRes[2]/M_PI*180. << "°" << std::endl;
	std::cout << "Sample delta two-theta = " << vecRes[1]/M_PI*180. << "°" << "\n";


	// generate plot
	std::ofstream ofstrGpl("powder.gpl");
	ofstrGpl << "bragg_2th(G, k, dtt) = 2.*asin(G/(2*k))/pi*180. + dtt\n";

	ofstrGpl << "\nset xlabel \"Wavenumber G (1/A)\"\n";
	ofstrGpl << "set ylabel \"Scattering angle (deg)\"\n";
	ofstrGpl << "unset key\n";

	ofstrGpl << "\nk = " << vecRes[0] << "\n";
	ofstrGpl << "dtt = " << vecRes[1]/M_PI*180. << "\n";

	ofstrGpl << "\nplot \\"
		<< "\n\t bragg_2th(x, k, dtt), \\"
		<< "\n\t\"-\" u 1:2:3 w yerrorbars pt 7\n";
	for(std::size_t i=0; i<vecGs.size(); ++i)
		ofstrGpl << vecGs[i] << "\t" << vecTTs[i]/M_PI*180. << "\t" << vecTTErrs[i]/M_PI*180. << "\n";
	ofstrGpl << "e\n";

	return 0;
}
