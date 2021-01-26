/**
 * tlibs test file
 * @author Tobias Weber <tobias.weber@tum.de>
 * @license GPLv2 or GPLv3
 */

// gcc -o coord coord.cpp -lstdc++ -std=c++11 -lm

#include "../math/math.h"
#include <iostream>
#include <fstream>

bool all_even(int h, int k, int l)
{
	return tl::is_even(h) && tl::is_even(k) && tl::is_even(l);
}

bool all_odd(int h, int k, int l)
{
	return tl::is_odd(h) && tl::is_odd(k) && tl::is_odd(l);
}

void to_gnomonic(double x, double y, double z, double& xg, double& yg)
{
	double rho, phi, theta, phi_crys, theta_crys;
	std::tie(rho,phi,theta) = tl::cart_to_sph(x, y, z);
	std::tie(phi_crys,theta_crys) = tl::sph_to_crys(phi, theta);
//	std::cout << "rho phi theta = " << rho << ", "
//		<< tl::r2d(phi) << ", " << tl::r2d(theta) << std::endl;
//	std::cout << "phi_cyrs theta_crys = "
//		<< tl::r2d(phi_crys) << ", " << tl::r2d(theta_crys) << std::endl;

//	std::tie(x,y,z) = tl::sph_to_cart(rho, phi, theta);
//	std::cout << "xyz = " << x << ", " << y << ", " << z << std::endl;

	std::tie(xg, yg) = tl::gnomonic_proj(phi_crys, theta_crys);
//	std::cout << "gnomonic xy = " << xg << ", " << yg << std::endl;
}


int main()
{
/*	while(1)
	{
		double phi, theta;
		std::cout << "phi theta = "; std::cin >> phi >> theta;
		phi = tl::d2r(phi);
		theta = tl::d2r(theta);

		double xg, yg;
		std::tie(xg, yg) = tl::gnomonic_proj(phi, theta);
		std::cout << "gnomonic xy = " << xg << ", " << yg << std::endl;
	}*/

	double drot = M_PI/4.;
	std::ofstream ofstr("/tmp/0.dat");

	for(int x=-7; x<=7; ++x)
	for(int y=-7; y<=7; ++y)
	for(int z=-7; z<=7; ++z)
	{
		if(!(all_even(x,y,z) || all_odd(x,y,z))) continue;

		double dx=x, dy=y, dz=z;

		dx = std::cos(drot)*dx - std::sin(drot)*dy;
		dy = std::sin(drot)*dx + std::cos(drot)*dy;

		double xg, yg;
		to_gnomonic(dx,dy,dz, xg,yg);

		ofstr << xg << "\t" << yg << "\t# " << x << "\t" << y << "\t" << z << "\n";
	}

	std::system("gnuplot -p -e \"set xrange[-1:1]; set yrange[-1:1]; plot '/tmp/0.dat' u 1:2\"");
	return 0;
}
