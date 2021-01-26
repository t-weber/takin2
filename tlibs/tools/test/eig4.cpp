/**
 * tlibs test file
 * @author Tobias Weber <tobias.weber@tum.de>
 * @license GPLv2 or GPLv3
 */

// gcc -I/usr/include/lapacke -o eig4 eig4.cpp ../math/linalg2.cpp ../log/log.cpp -lstdc++ -lm -llapacke -llapack -std=c++11

#define TLIBS_INC_HDR_IMPLS
#include "../math/linalg.h"
#include "../math/linalg2.h"

namespace ublas = boost::numeric::ublas;
using T = double;

int main()
{
	ublas::matrix<T> M = tl::make_mat({
		{1.1,2.2,3.0},
		{2.2,5.5,6.6},
		{3.3,6.6,9.9}	});

	ublas::matrix<T> M_org = M;
	std::cout << "M =" << M << std::endl;
	std::cout << "symm: " << tl::is_symmetric(M) << std::endl;

	std::vector<ublas::vector<T>> evecs2_r, evecs2_i;
	std::vector<T> evals2_r, evals2_i;
	tl::eigenvec(M, evecs2_r, evecs2_i, evals2_r, evals2_i);
	for(int i=0; i<evals2_r.size(); ++i)
		std::cout << "eval r: " << evals2_r[i] <<
		", evec r: " << evecs2_r[i] << std::endl;
	for(int i=0; i<evals2_i.size(); ++i)
		std::cout << "eval i: " << evals2_i[i] <<
		", evec i: " << evecs2_i[i] << std::endl;
	std::cout << std::endl;


	ublas::matrix<T> MtM = ublas::prod(ublas::trans(M), M);
	std::cout << "M^tM = " << MtM << std::endl;
	std::cout << "symm: " << tl::is_symmetric(MtM) << std::endl;
	std::vector<ublas::vector<T>> evecs3_r, evecs3_i;
	std::vector<T> evals3_r, evals3_i;
	tl::eigenvec(MtM, evecs3_r, evecs3_i, evals3_r, evals3_i);
	for(int i=0; i<evals3_r.size(); ++i)
		std::cout << "eval r: " << evals3_r[i] <<
		", sqrt: " << std::sqrt(std::abs(evals3_r[i])) <<
		", evec r: " << evecs3_r[i] << std::endl;
	for(int i=0; i<evals3_i.size(); ++i)
		std::cout << "eval i: " << evals3_i[i] <<
		", evec i: " << evecs3_i[i] << std::endl;
	std::cout << std::endl;


	std::vector<ublas::vector<T>> evecs4;;
	std::vector<T> evals4;
	tl::eigenvec_approxsym/*_simple*/(M, evecs4, evals4);
	for(int i=0; i<evals4.size(); ++i)
		std::cout << "eval: " << evals4[i] <<
		", evec: " << evecs4[i] << std::endl;

	return 0;
}
