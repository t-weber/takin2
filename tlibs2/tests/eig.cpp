/**
 * test eigensystem calculation
 * @author Tobias Weber <tweber@ill.fr>
 * @date 25-jul-20
 * @license GPLv3, see 'LICENSE' file
 *
 * g++ -std=c++20 -DUSE_LAPACK -I.. -I/usr/include/lapacke -Iext/lapacke/include -Lext/lapacke/lib -o eig eig.cpp -llapacke
 */

#define BOOST_TEST_MODULE Eigenvector Test
#include <boost/test/included/unit_test.hpp>
namespace test = boost::unit_test;
namespace testtools = boost::test_tools;

#include <iostream>
#include <vector>
#include <random>

#include "libs/math20.h"


// LinearAlgebra.eigen([-1.5 0.01 0.02; -0.04 1.0 0.03; -2.05 0.06 0.5])

using t_types = std::tuple<double, float>;
BOOST_AUTO_TEST_CASE_TEMPLATE(test_eig, t_real, t_types)
{
	using namespace tl2_ops;

	using t_cplx = std::complex<t_real>;
	using t_vec = tl2::vec<t_real, std::vector>;
	using t_mat = tl2::mat<t_real, std::vector>;
	using t_vec_cplx = tl2::vec<t_cplx, std::vector>;
	using t_mat_cplx = tl2::mat<t_cplx, std::vector>;

	t_real eps = 1e-4;
	std::size_t dim = 3;

	std::mt19937 rndgen{tl2::epoch<unsigned int>()};
	std::uniform_real_distribution<t_real> rnddist{-100, 100};

	// real version
	{
		std::cout << "--------------------------------------------------------------------------------" << std::endl;
		auto mat = tl2::zero<t_mat>(dim,dim);
		for(std::size_t i=0; i<dim; ++i)
			for(std::size_t j=0; j<dim; ++j)
				mat(i,j) = rnddist(rndgen);
		std::cout << mat << std::endl;

		bool sym = 0;
		auto [ok, evals_re, evals_im, evecs_re, evecs_im] =
			tl2_la::eigenvec<t_mat, t_vec, t_real>(mat, false, sym, true);
		std::cout << "ok = " << std::boolalpha << ok << "\n" << std::endl;

		for(std::size_t i=0; i<evals_re.size(); ++i)
			std::cout << "Re(eval) " << i+1 << ": " << evals_re[i] << std::endl;
		for(std::size_t i=0; i<evals_im.size(); ++i)
			std::cout << "Im(eval) " << i+1 << ": " << evals_im[i] << std::endl;
		std::cout << std::endl;

		for(std::size_t i=0; i<evecs_re.size(); ++i)
			std::cout << "Re(evec) " << i+1 << ": " << evecs_re[i] << std::endl;
		std::cout << std::endl;
		for(std::size_t i=0; i<evecs_im.size(); ++i)
			std::cout << "Im(evec) " << i+1 << ": " << evecs_im[i] << std::endl;
		std::cout << std::endl;

		BOOST_TEST(ok);


		t_mat_cplx mat_cplx = tl2::zero<t_mat_cplx>(dim, dim);
		for(std::size_t i=0; i<dim; ++i)
			for(std::size_t j=0; j<dim; ++j)
				mat_cplx(i,j) = mat(i,j);

		for(std::size_t i=0; i<dim; ++i)
		{
			t_cplx eval = t_cplx{evals_re[i], evals_im[i]};
			t_vec_cplx evec = tl2::zero<t_vec_cplx>(dim);

			for(std::size_t j=0; j<dim; ++j)
				evec[j] = t_cplx{evecs_re[i][j], evecs_im[i][j]};

			auto tstvec1 = mat_cplx * evec;
			auto tstvec2 = eval * evec;
			bool is_equal = tl2::equals<t_vec_cplx>(tstvec1, tstvec2, eps);
			std::cout << tstvec1 << " == " << tstvec2 << ": " << std::boolalpha << is_equal << std::endl;
			BOOST_TEST(is_equal);
		}

		std::cout << "--------------------------------------------------------------------------------\n" << std::endl;
	}


	// real version, symmetric
	{
		std::cout << "--------------------------------------------------------------------------------" << std::endl;
		auto mat = tl2::zero<t_mat>(dim,dim);
		for(std::size_t i=0; i<dim; ++i)
			for(std::size_t j=i; j<dim; ++j)
				mat(j,i) = mat(i,j) = rnddist(rndgen);
		std::cout << mat << std::endl;

		bool sym = 1;
		auto [ok, evals_re, evals_im, evecs_re, evecs_im] =
			tl2_la::eigenvec<t_mat, t_vec, t_real>(mat, false, sym, true);
		std::cout << "ok = " << std::boolalpha << ok << "\n" << std::endl;

		for(std::size_t i=0; i<evals_re.size(); ++i)
			std::cout << "Re(eval) " << i+1 << ": " << evals_re[i] << std::endl;
		for(std::size_t i=0; i<evals_im.size(); ++i)
			std::cout << "Im(eval) " << i+1 << ": " << evals_im[i] << std::endl;
		std::cout << std::endl;

		for(std::size_t i=0; i<evecs_re.size(); ++i)
			std::cout << "Re(evec) " << i+1 << ": " << evecs_re[i] << std::endl;
		std::cout << std::endl;
		for(std::size_t i=0; i<evecs_im.size(); ++i)
			std::cout << "Im(evec) " << i+1 << ": " << evecs_im[i] << std::endl;
		std::cout << std::endl;

		BOOST_TEST(ok);


		t_mat_cplx mat_cplx = tl2::zero<t_mat_cplx>(dim, dim);
		for(std::size_t i=0; i<dim; ++i)
			for(std::size_t j=0; j<dim; ++j)
				mat_cplx(i,j) = mat(i,j);

		for(std::size_t i=0; i<dim; ++i)
		{
			t_cplx eval = t_cplx{evals_re[i], evals_im[i]};
			t_vec_cplx evec = tl2::zero<t_vec_cplx>(dim);

			for(std::size_t j=0; j<dim; ++j)
				evec[j] = t_cplx{evecs_re[i][j], evecs_im[i][j]};

			auto tstvec1 = mat_cplx * evec;
			auto tstvec2 = eval * evec;
			bool is_equal = tl2::equals<t_vec_cplx>(tstvec1, tstvec2, eps);
			std::cout << tstvec1 << " == " << tstvec2 << ": " << std::boolalpha << is_equal << std::endl;
			BOOST_TEST(is_equal);
		}

		std::cout << "--------------------------------------------------------------------------------\n" << std::endl;
	}


	// complex version
	{
		std::cout << "--------------------------------------------------------------------------------" << std::endl;

		auto mat = tl2::zero<t_mat_cplx>(dim, dim);
		for(std::size_t i=0; i<dim; ++i)
			for(std::size_t j=0; j<dim; ++j)
				mat(i,j) = rnddist(rndgen) + rnddist(rndgen)*t_cplx{0, 1};
		std::cout << mat << std::endl;

		bool herm = 0;
		auto [ok, evals, evecs] =
			tl2_la::eigenvec<t_mat_cplx, t_vec_cplx, t_cplx>(mat, false, herm, true);
		std::cout << "ok = " << std::boolalpha << ok << "\n" << std::endl;

		for(std::size_t i=0; i<evecs.size(); ++i)
			std::cout << "eval " << i+1 << ": " << evals[i] << std::endl;
		std::cout << std::endl;
		for(std::size_t i=0; i<evecs.size(); ++i)
			std::cout << "evec " << i+1 << ": " << evecs[i] << std::endl;
		std::cout << std::endl;

		BOOST_TEST(ok);
		for(std::size_t i=0; i<dim; ++i)
		{
			auto tstvec1 = mat*evecs[i];
			auto tstvec2 = evals[i]*evecs[i];
			bool is_equal = tl2::equals<t_vec_cplx>(tstvec1, tstvec2, eps);
			std::cout << tstvec1 << " == " << tstvec2 << ": " << std::boolalpha << is_equal << std::endl;
			BOOST_TEST(is_equal);
		}

		std::cout << "--------------------------------------------------------------------------------\n" << std::endl;
	}


	// complex version, hermitian
	{
		std::cout << "--------------------------------------------------------------------------------" << std::endl;

		auto mat = tl2::zero<t_mat_cplx>(dim, dim);
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

		bool herm = 1;
		auto [ok, evals, evecs] =
			tl2_la::eigenvec<t_mat_cplx, t_vec_cplx, t_cplx>(mat, false, herm, true);
		std::cout << "ok = " << std::boolalpha << ok << "\n" << std::endl;

		for(std::size_t i=0; i<evecs.size(); ++i)
			std::cout << "eval " << i+1 << ": " << evals[i] << std::endl;
		std::cout << std::endl;
		for(std::size_t i=0; i<evecs.size(); ++i)
			std::cout << "evec " << i+1 << ": " << evecs[i] << std::endl;
		std::cout << std::endl;

		BOOST_TEST(ok);
		for(std::size_t i=0; i<dim; ++i)
		{
			auto tstvec1 = mat*evecs[i];
			auto tstvec2 = evals[i]*evecs[i];
			bool is_equal = tl2::equals<t_vec_cplx>(tstvec1, tstvec2, eps);
			std::cout << tstvec1 << " == " << tstvec2 << ": " << std::boolalpha << is_equal << std::endl;
			BOOST_TEST(is_equal);
		}

		std::cout << "--------------------------------------------------------------------------------\n" << std::endl;
	}
}
