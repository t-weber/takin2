/**
 * advanced linalg helpers
 * @author: Tobias Weber <tobias.weber@tum.de>
 * @date: 30-apr-2013 - 2018
 * @license GPLv2 or GPLv3
 */

#ifndef __TLIBS_LINALG2__
#define __TLIBS_LINALG2__


#include "math.h"
#include "linalg.h"
#include <complex>


namespace tl {

#if !defined NO_LAPACK && !defined USE_LAPACK
	#define USE_LAPACK
#endif


#ifdef NO_LAPACK	// direct implementation

/**
 * qr decomposition: M = QR
 */
template<typename T/*=double*/>
bool qr(const ublas::matrix<T>& M, ublas::matrix<T>& Q, ublas::matrix<T>& R)
{
	return qr_decomp(M, Q, R);
}

/**
 * calculates the eigenvectors of a symmetric matrix
 */
template<typename T=double>
bool eigenvec_sym(const ublas::matrix<T>& mat,
	std::vector<ublas::vector<T>>& evecs, std::vector<T>& evals, bool bNorm=0)
{
	bool bOk = eigenvec_sym_simple(mat, evecs, evals);

	if(bNorm && bOk)
	{
		for(std::size_t i=0; i<evecs.size(); ++i)
			evecs[i] /= veclen(evecs[i]);
	}

	return bOk;
}

/**
 * calculates the approximate eigenvectors
 */
template<typename T=double>
bool eigenvec_approxsym(const ublas::matrix<T>& mat,
	std::vector<ublas::vector<T>>& evecs, std::vector<T>& evals, bool bNorm=0)
{
	bool bOk = eigenvec_approxsym_simple(mat, evecs, evals);

	if(bNorm && bOk)
	{
		for(std::size_t i=0; i<evecs.size(); ++i)
			evecs[i] /= veclen(evecs[i]);
	}

	return bOk;
}


#else	// Lapack wrappers


/**
 * qr decomposition: M = QR
 */
template<typename T/*=double*/>
bool qr(const ublas::matrix<T>& M, ublas::matrix<T>& Q, ublas::matrix<T>& R);


/**
 * calculates the eigenvectors of a general matrix
 */
template<typename T=double>
bool eigenvec(const ublas::matrix<T>& mat,
	std::vector<ublas::vector<T> >& evecs_real, std::vector<ublas::vector<T>>& evecs_imag,
	std::vector<T>& evals_real, std::vector<T>& evals_imag, bool bNorm=0);

/**
 * calculates only the eigenvalues of a general matrix
 */
template<typename T=double>
bool eigenval(const ublas::matrix<T>& mat,
	std::vector<T>& evals_real, std::vector<T>& evals_imag);

/**
 * calculates the eigenvectors of a general complex matrix
 */
template<typename T=double>
bool eigenvec_cplx(const ublas::matrix<std::complex<T>>& mat,
	std::vector<ublas::vector<std::complex<T>> >& evecs,
	std::vector<std::complex<T>>& evals, bool bNorm=0);

/**
 * calculates only the eigenvalues of a general complex matrix
 */
template<typename T=double>
bool eigenval_cplx(const ublas::matrix<std::complex<T>>& mat, std::vector<std::complex<T>>& evals);


/**
 * calculates the eigenvectors of a symmetric matrix
 */
template<typename T=double>
bool eigenvec_sym(const ublas::matrix<T>& mat,
	std::vector<ublas::vector<T>>& evecs, std::vector<T>& evals, bool bNorm=0);

/**
 * calculates only the eigenvalues of a symmetric matrix
 */
template<typename T=double>
bool eigenval_sym(const ublas::matrix<T>& mat, std::vector<T>& evals);

/**
 * calculates the eigenvectors of a hermitian matrix
 */
template<typename T=double>
bool eigenvec_herm(const ublas::matrix<std::complex<T>>& mat,
	std::vector<ublas::vector<std::complex<T>>>& evecs,
	std::vector<T>& evals, bool bNorm=0);

/**
 * calculates only the eigenvalues of a hermitian matrix
 */
template<typename T=double>
bool eigenval_herm(const ublas::matrix<std::complex<T>>& mat, std::vector<T>& evals);


/**
 * calculates selected eigenvectors of a hermitian matrix
 */
template<typename T=double>
bool eigenvecsel_herm(const ublas::matrix<std::complex<T>>& mat,
	std::vector<ublas::vector<std::complex<T>>>& evecs,
	std::vector<T>& evals, bool bNorm=0, T minval=-1, T maxval=-2, T eps=T(-1));


/**
 * calculates the singular values of a real matrix: M = U diag(vals) V^t
 */
template<typename T=double>
bool singvec(const ublas::matrix<T>& mat,
	ublas::matrix<T>& matU, ublas::matrix<T>& matV, std::vector<T>& vecsvals);

/**
 * calculates the singular values of a complex matrix: M = U diag(vals) (V*)^t
 */
template<typename T=double>
bool singvec_cplx(const ublas::matrix<std::complex<T>>& mat,
	ublas::matrix<std::complex<T>>& matU, ublas::matrix<std::complex<T>>& matV,
	std::vector<T>& vecsvals);

/**
 * pseudoinverse of a real diagonal matrix
 * see: https://de.wikipedia.org/wiki/Pseudoinverse#Berechnung
 */
template<typename T=double>
ublas::matrix<T> pseudoinverse_diag(const ublas::matrix<T>& mat)
{
	std::size_t N = std::min(mat.size1(), mat.size2());
	ublas::matrix<T> matRet = mat;

	for(std::size_t i=0; i<N; ++i)
	{
		if(!float_equal(mat(i,i), T(0)))
			matRet(i,i) = T(1) / mat(i,i);
	}

	return matRet;
}

/**
 * pseudoinverse M+ of a real matrix
 * M  = U D (V*)^t
 * M+ = V D+ (U*)^t
 * see: https://de.wikipedia.org/wiki/Pseudoinverse#Berechnung
 */
template<typename T=double>
bool pseudoinverse(const ublas::matrix<T>& mat, ublas::matrix<T>& matInv)
{
	std::vector<T> vecS;
	ublas::matrix<T> matU, matV;

	if(!singvec(mat, matU, matV, vecS))
		return false;

	ublas::matrix<T> matS = diag_matrix(vecS);
	matS = pseudoinverse_diag(matS);
	matU = transpose(matU);

	matInv = prod_mm(matS, matU);
	matInv = prod_mm(matV, matInv);

	return true;
}


/**
 * calculates the approximate eigenvectors
 */
template<typename T=double>
bool eigenvec_approxsym(const ublas::matrix<T>& mat,
	std::vector<ublas::vector<T>>& evecs, std::vector<T>& evals)
{
	ublas::matrix<T> matU, matV;
	bool bOk = singvec(mat, matU, matV, evals);

	evecs.resize(matV.size2());
	for(std::size_t j=0; j<matV.size2(); ++j)
		evecs[j] = get_column(matV, j);

	return bOk;
}



#ifdef TLIBS_INC_HDR_IMPLS
}
#include "linalg2_impl.h"
namespace tl {
#endif


#endif

}

#endif
