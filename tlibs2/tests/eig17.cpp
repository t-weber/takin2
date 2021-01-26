/**
 * test eigensystem calculation
 * @author Tobias Weber <tweber@ill.fr>
 * @date 11-aug-20
 * @license GPLv3, see 'LICENSE' file
 *
 * g++ -std=c++17 -DUSE_LAPACK -I.. -I/usr/include/lapacke -Iext/lapacke/include -Lext/lapacke/lib -o eig17 eig17.cpp ../libs/log.cpp -llapacke
 */

#include <iostream>
#include <vector>
#include <random>

#include "libs/math17.h"


using t_real = double;
using t_mat = ublas::matrix<t_real>;
using t_vec = ublas::vector<t_real>;
using t_cplx = std::complex<double>;
using t_mat_cplx = ublas::matrix<t_cplx>;
using t_vec_cplx = ublas::vector<t_cplx>;


std::size_t tst_real()
{
	t_real eps = 1e-4;
	std::size_t dim = 3;

	std::mt19937 rndgen{tl2::epoch<unsigned int>()};
	std::uniform_real_distribution<t_real> rnddist{-100, 100};

	auto mat = tl2::zero_matrix<t_mat>(dim, dim);
	for(std::size_t i=0; i<dim; ++i)
		for(std::size_t j=0; j<dim; ++j)
			mat(i,j) = rnddist(rndgen);

	std::cout << mat << std::endl;


	std::vector<t_vec> evecs_re, evecs_im;
	std::vector<t_real> evals_re, evals_im;

	std::cout << std::boolalpha << tl2::eigenvec(
		mat, evecs_re, evecs_im, evals_re, evals_im, false)
		<< std::endl;

	for(t_real val : evals_re)
		std::cout << "Re{eval}: " << val << std::endl;
	for(t_real val : evals_im)
		std::cout << "Im{eval}: " << val << std::endl;
	for(const t_vec& vec : evecs_re)
		std::cout << "Re{evec}: " << vec << std::endl;
	for(const t_vec& vec : evecs_im)
		std::cout << "Im{evec}: " << vec << std::endl;

	t_mat_cplx mat_cplx = tl2::zero_matrix<t_mat_cplx>(dim, dim);
	for(std::size_t i=0; i<dim; ++i)
		for(std::size_t j=0; j<dim; ++j)
			mat_cplx(i,j) = mat(i,j);

	std::size_t failures = 0;
	for(std::size_t i=0; i<dim; ++i)
	{
		t_cplx eval = t_cplx{evals_re[i], evals_im[i]};
		t_vec_cplx evec = tl2::zero_vector<t_vec_cplx>(dim);

		for(std::size_t j=0; j<dim; ++j)
			evec[j] = t_cplx{evecs_re[i][j], evecs_im[i][j]};

		auto tstvec1 = tl2::prod_mv(mat_cplx, evec);
		auto tstvec2 = eval * evec;
		bool is_equal = tl2::vec_equal<t_vec_cplx>(tstvec1, tstvec2, eps);
		std::cout << tstvec1 << " == " << tstvec2 << ": " << std::boolalpha << is_equal << std::endl;
	}

	return failures;
}


std::size_t tst_cplx()
{
	t_real eps = 1e-4;
	std::size_t dim = 3;

	std::mt19937 rndgen{tl2::epoch<unsigned int>()};
	std::uniform_real_distribution<t_real> rnddist{-100, 100};

	auto mat = tl2::zero_matrix<t_mat_cplx>(dim, dim);
	for(std::size_t i=0; i<dim; ++i)
		for(std::size_t j=0; j<dim; ++j)
			mat(i,j) = rnddist(rndgen) + rnddist(rndgen)*t_cplx{0, 1};

	std::cout << mat << std::endl;


	std::vector<t_vec_cplx> evecs;
	std::vector<t_cplx> evals;

	std::cout << std::boolalpha << tl2::eigenvec_cplx(mat, evecs, evals, false) << std::endl;

	for(const t_cplx& val : evals)
		std::cout << "eval: " << val << std::endl;
	for(const t_vec_cplx& vec : evecs)
		std::cout << "evec: " << vec << std::endl;

	std::size_t failures = 0;
	for(std::size_t i=0; i<dim; ++i)
	{
		auto tstvec1 = tl2::prod_mv(mat, evecs[i]);
		auto tstvec2 = evals[i]*evecs[i];
		bool is_equal = tl2::vec_equal<t_vec_cplx>(tstvec1, tstvec2, eps);
		if(!is_equal)
			++failures;
		std::cout << tstvec1 << " == " << tstvec2 << ": " << std::boolalpha << is_equal << std::endl;
	}

	return failures;
}


std::size_t tst_herm()
{
	t_real eps = 1e-4;
	std::size_t dim = 3;

	std::mt19937 rndgen{tl2::epoch<unsigned int>()};
	std::uniform_real_distribution<t_real> rnddist{-100, 100};

	auto mat = tl2::zero_matrix<t_mat_cplx>(dim, dim);
	for(std::size_t i=0; i<dim; ++i)
	{
		for(std::size_t j=i+1; j<dim; ++j)
		{
			mat(i,j) = rnddist(rndgen) + rnddist(rndgen)*t_cplx{0, 1};
			mat(j,i) = std::conj(mat(i,j));
		}
	}
	for(std::size_t i=0; i<dim; ++i)
		mat(i,i) = rnddist(rndgen);

	std::cout << mat << std::endl;


	std::vector<t_vec_cplx> evecs;
	std::vector<t_real> evals;

	std::cout << std::boolalpha << tl2::eigenvec_herm(mat, evecs, evals, false) << std::endl;

	for(t_real val : evals)
		std::cout << "eval: " << val << std::endl;
	for(const t_vec_cplx& vec : evecs)
		std::cout << "evec: " << vec << std::endl;

	std::size_t failures = 0;
	for(std::size_t i=0; i<dim; ++i)
	{
		auto tstvec1 = tl2::prod_mv(mat, evecs[i]);
		auto tstvec2 = evals[i]*evecs[i];
		bool is_equal = tl2::vec_equal<t_vec_cplx>(tstvec1, tstvec2, eps);
		if(!is_equal)
			++failures;
		std::cout << tstvec1 << " == " << tstvec2 << ": " << std::boolalpha << is_equal << std::endl;
	}

	return failures;
}


std::size_t tst_herm_sel()
{
	t_real eps = 1e-4;
	std::size_t dim = 3;

	std::mt19937 rndgen{tl2::epoch<unsigned int>()};
	std::uniform_real_distribution<t_real> rnddist{-100, 100};

	auto mat = tl2::zero_matrix<t_mat_cplx>(dim, dim);
	for(std::size_t i=0; i<dim; ++i)
	{
		for(std::size_t j=i+1; j<dim; ++j)
		{
			mat(i,j) = rnddist(rndgen) + rnddist(rndgen)*t_cplx{0, 1};
			mat(j,i) = std::conj(mat(i,j));
		}
	}
	for(std::size_t i=0; i<dim; ++i)
		mat(i,i) = rnddist(rndgen);

	std::cout << mat << std::endl;


	std::vector<t_vec_cplx> evecs;
	std::vector<t_real> evals;

	std::cout << std::boolalpha << tl2::eigenvecsel_herm(mat, evecs, evals, false) << std::endl;

	for(t_real val : evals)
		std::cout << "eval: " << val << std::endl;
	for(const t_vec_cplx& vec : evecs)
		std::cout << "evec: " << vec << std::endl;

	std::size_t failures = 0;
	for(std::size_t i=0; i<dim; ++i)
	{
		auto tstvec1 = tl2::prod_mv(mat, evecs[i]);
		auto tstvec2 = evals[i]*evecs[i];
		bool is_equal = tl2::vec_equal<t_vec_cplx>(tstvec1, tstvec2, eps);
		if(!is_equal)
			++failures;
		std::cout << tstvec1 << " == " << tstvec2 << ": " << std::boolalpha << is_equal << std::endl;
	}

	return failures;
}


int main()
{
	std::cout << "General real matrix" << std::endl;
	std::size_t failures = tst_real();
	std::cerr << "Failures: " << failures << std::endl;

	std::cout << "\nGeneral complex matrix" << std::endl;
	failures = tst_cplx();
	std::cerr << "Failures: " << failures << std::endl;

	std::cout << "\nHermitian, selected values" << std::endl;
	failures = tst_herm_sel();
	std::cerr << "Failures: " << failures << std::endl;

	// there seems to be a bug in lapacke: https://github.com/Reference-LAPACK/lapack/issues/379
	std::cout << "\nHermitian" << std::endl;
	failures = tst_herm();
	std::cerr << "Failures: " << failures << std::endl;

	return 0;
}
