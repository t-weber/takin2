/**
 * tlibs test file
 * @author Tobias Weber <tobias.weber@tum.de>
 * @license GPLv2 or GPLv3
 */

// gcc -I/usr/include/lapacke -o eig5 eig5.cpp ../math/linalg2.cpp ../log/log.cpp -lstdc++ -lm -llapacke -llapack -std=c++11

#define TLIBS_INC_HDR_IMPLS
#include "../math/linalg.h"
#include "../math/linalg2.h"

namespace ublas = boost::numeric::ublas;
using T = double;

template
bool tl::singvec_cplx(const ublas::matrix<std::complex<double>>& mat,
    ublas::matrix<std::complex<double>>& matU, ublas::matrix<std::complex<double>>& matV,
    std::vector<double>& vecsvals);
template
bool tl::singvec_cplx(const ublas::matrix<std::complex<float>>& mat,
    ublas::matrix<std::complex<float>>& matU, ublas::matrix<std::complex<float>>& matV,
    std::vector<float>& vecsvals);


int main()
{
	ublas::matrix<T> M = tl::make_mat({
		{1.1,2.2,3.0},
		{2.2,5.5,6.6},
		{3.3,6.6,9.9}	});

	ublas::matrix<T> M_org = M;
	std::cout << "M = " << M << std::endl;
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


	std::vector<ublas::vector<T>> evecs4;
	std::vector<T> evals4;
	tl::eigenvec_approxsym_simple(M, evecs4, evals4);
	for(int i=0; i<evals4.size(); ++i)
		std::cout << "eval: " << evals4[i] <<
		", evec: " << evecs4[i] << std::endl;
	std::cout << std::endl;


	std::vector<ublas::vector<T>> evecs5;
	std::vector<T> evals5;
	tl::eigenvec_approxsym(M, evecs5, evals5);
	for(int i=0; i<evals5.size(); ++i)
		std::cout << "eval: " << evals5[i] <<
		", evec: " << evecs5[i] << std::endl;
	std::cout << std::endl;


	ublas::matrix<T> matU, matV;
	std::vector<T> vecSvals;
	tl::singvec(M, matU, matV, vecSvals);
	std::cout << "U = " << matU << std::endl;
	std::cout << "V = " << matV << std::endl;
	std::cout << "svals: ";
	for(int i=0; i<vecSvals.size(); ++i)
		std::cout << vecSvals[i] << ", ";
	std::cout << std::endl;

	ublas::matrix<T> matD = tl::diag_matrix(vecSvals);
	ublas::matrix<T> matDVt = ublas::prod(matD, ublas::trans(matV));
	std::cout << "U*D*V^t = " << ublas::prod(matU, matDVt) << std::endl;
	std::cout << std::endl;

	ublas::matrix<T> Mp;
	tl::pseudoinverse(M, Mp);
	std::cout << "M+ = " << Mp << std::endl;

	return 0;
}
