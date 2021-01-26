/**
 * tlibs test file
 * @author Tobias Weber <tobias.weber@tum.de>
 * @license GPLv2 or GPLv3
 */

// clang -o hund hund.cpp -std=c++11 -lstdc++

#include <iostream>
#include "../math/term.h"

int main()
{
	//auto tup = tl::get_orbital<std::string>("2p3");
	//std::cout << std::get<0>(tup) << ", " << std::get<1>(tup) << ", " << std::get<2>(tup) << std::endl;


	double S,L,J;
	std::uint16_t l, es;

	std::cout << "l = "; std::cin >> l;
	std::cout << "Number of electrons: "; std::cin >> es;

	//std::tie(S,L,J) = tl::hund(std::string("1s2 2p3"));
	std::tie(S,L,J) = tl::hund(l, es);
	std::string strTerm = tl::get_termsymbol(S,L,J);

	std::cout << strTerm << "\n";
	std::cout << "S = " << S << "\n";
	std::cout << "L = " << L << "\n";
	std::cout << "J = " << J << "\n";
	std::cout.flush();

	return 0;
}
