/**
 * tensor function
 * @author Tobias Weber <tobias.weber@tum.de>
 * @date may-2017
 * @license GPLv2 or GPLv3
 */

#ifndef __TLIBS_TENSOR_H__
#define __TLIBS_TENSOR_H__


#include "math.h"
#include "linalg.h"
#include "linalg_ops.h"
#include <algorithm>
#include <numeric>


namespace tl {


/**
 * count permutations to get from "idx1" to "idx2" (idx2 has to be sorted)
 * WARNING: doesn't work: next_permutation randomly shuffles the permutations!
 */
template<class t_lst = std::vector<std::size_t>>
std::pair<bool, std::size_t> count_permutations(const t_lst& idx1, const t_lst& _idx2)
{
	std::size_t iPerms = 0;
	bool bHasPerm = 1;
	t_lst idx2 = _idx2;

	while(1)
	{
		if(std::equal(idx1.begin(), idx1.end(), idx2.begin()))
			break;
		if(!std::next_permutation(idx2.begin(), idx2.end()))
		{
			bHasPerm = 0;
			break;
		}

		++iPerms;
	}

	return std::make_pair(bHasPerm, iPerms);
}


/**
 * count permutations to get to ordered indices (0-based)
 */
template<class t_lst = std::vector<std::size_t>>
std::pair<bool, std::size_t> count_permutations(const t_lst& _vecIdx)
{
	std::size_t iNumPerms = 0;
	bool bHasPerm = 1;
	t_lst vecIdx = _vecIdx;

	for(std::size_t iIdx = 0; iIdx < vecIdx.size(); ++iIdx)
	{
		// already at correct position?
		if(vecIdx[iIdx] == iIdx)
			continue;

		// else permutate
		auto iterPos = std::find(vecIdx.begin()+iIdx, vecIdx.end(), iIdx);
		if(iterPos == vecIdx.end())
		{
			// no permutation found
			bHasPerm = 0;
			break;
		}

		std::iter_swap(vecIdx.begin()+iIdx, iterPos);
		++iNumPerms;
	}

	return std::make_pair(bHasPerm, iNumPerms);
}


/**
 * elements of the (cartesian) epsilon tensor (indices 0-based)
 */
template<typename T = double, class t_lst = std::vector<std::size_t>>
T epsilon_tensor(const t_lst& idx)
{
	//t_lst idxIota(idx.size());
	//std::size_t iMin = *std::min_element(idx.begin(), idx.end());
	//std::iota(idxIota.begin(), idxIota.end(), iMin);

	std::size_t iNumPerms;
	bool bHasPerm;
	std::tie(bHasPerm, iNumPerms) = count_permutations<t_lst>(idx/*, idxIota*/);

	if(!bHasPerm)
		return T(0);
	if(is_even(iNumPerms))
		return T(1);
	return T(-1);
}


/**
 * elements of the (non-cartesian) epsilon tensor (indices 0-based)
 */
template<typename T = double, class t_lst = std::vector<std::size_t>,
	class t_mat = ublas::matrix<double>>
T epsilon_tensor(const t_mat& matGcov, const t_lst& idx, bool bCov = 1)
{
	T tEps = epsilon_tensor<T, t_lst>(idx);
	T tMetr = std::abs(determinant(matGcov));
	T tDet = std::sqrt(tMetr)*tEps;

	if(!bCov && !float_equal(tDet, T(0)))
		tDet = T(1) / tDet;
	return tDet;
}


/**
 * creates a metric tensor
 */
template<class t_mat = ublas::matrix<double>,
	class t_vec = ublas::vector<double>,
	template<class...> class t_lst = std::initializer_list>
t_mat make_metric_cov(const t_lst<t_vec>& lstVecsCov)
{
	using T = typename t_mat::value_type;
	const std::size_t iDim = std::min(lstVecsCov.size(), lstVecsCov.begin()->size());
	t_mat matG(iDim, iDim);

	typename t_lst<t_vec>::const_iterator iter;
	std::size_t i;
	for(iter=lstVecsCov.begin(), i=0; i<iDim; ++iter, ++i)
	{
		typename t_lst<t_vec>::const_iterator iter2;
		std::size_t j;
		for(iter2=lstVecsCov.begin(), j=0; j<iDim; ++iter2, ++j)
		{
			matG(i,j) = mult<t_vec, t_vec>(*iter, *iter2);
		}
	}

	return matG;
}


/**
 * inner product using metric matGCov
 */
template<class t_mat = ublas::matrix<double>,
	class t_vec = ublas::vector<double>>
typename t_vec::value_type inner_prod(const t_mat& matGCov,
	const t_vec& vec1Contra, const t_vec& vec2Contra)
{
	t_vec vec1Cov = mult<t_mat, t_vec>(matGCov, vec1Contra);
	return mult<t_vec, t_vec>(vec1Cov, vec2Contra);
}


/**
 * vector length
 */
template<class t_mat = ublas::matrix<double>,
	class t_vec = ublas::vector<double>>
typename t_vec::value_type vec_len(const t_mat& matGCov,
	const t_vec& vecContra)
{
	using T = typename t_vec::value_type;

	T tdot = inner_prod(matGCov, vecContra, vecContra);
	return std::sqrt(tdot);
}


/**
 * angle between vectors
 */
template<class t_mat = ublas::matrix<double>,
	class t_vec = ublas::vector<double>>
typename t_vec::value_type vec_angle(const t_mat& matGCov,
	const t_vec& vec1Contra, const t_vec& vec2Contra)
{
	using T = typename t_vec::value_type;

	T lenv1 = vec_len<t_mat, t_vec>(matGCov, vec1Contra);
	T lenv2 = vec_len<t_mat, t_vec>(matGCov, vec2Contra);
	T v1v2 = inner_prod<t_mat, t_vec>(matGCov, vec1Contra, vec2Contra);

	return std::acos(v1v2 / (lenv1 * lenv2));
}


/**
 * cross product using metric matGContra
 * see: (Arens 2015), p. 815
 */
template<class t_mat = ublas::matrix<double>,
	class t_vec = ublas::vector<double>>
t_vec cross_prod_contra(const t_mat& matGCov,
	const t_vec& vec1Contra, const t_vec& vec2Contra,
	bool bNorm = false)
{
	using T = typename t_mat::value_type;
	const std::size_t iDim = matGCov.size1();

	t_mat matGContra;
	inverse(matGCov, matGContra);

	t_vec vecCrossContra = zero_v<t_vec>(iDim);

	for(std::size_t j=0; j<iDim; ++j)
	{
		for(std::size_t k=0; k<iDim; ++k)
		{
			for(std::size_t l=0; l<iDim; ++l)
			{
				T eps = epsilon_tensor(matGCov, {l,j,k}, true);
				if(eps == T(0)) continue;

				for(std::size_t i=0; i<iDim; ++i)
					vecCrossContra[i] += matGContra(i,l)*eps * vec1Contra[j]*vec2Contra[k];
			}
		}
	}

	if(bNorm)
		vecCrossContra /= vec_len<t_mat, t_vec>(matGCov, vecCrossContra);
	return vecCrossContra;
}



/**
 * tensor product
 * see e.g.: (Arfken 2013), p. 109
 */
template<class t_mat = ublas::matrix<double>>
t_mat tensor_prod(const t_mat& mat1, const t_mat& mat2)
{
	t_mat mat(mat1.size1()*mat2.size1(), mat1.size2()*mat2.size2());

	for(std::size_t i=0; i<mat1.size1(); ++i)
	{
		for(std::size_t j=0; j<mat1.size2(); ++j)
		{
			t_mat matElem = mat1(i,j)*mat2;
			submatrix_copy(mat, matElem,
				i*mat2.size1(), j*mat2.size2());
		}
	}
	return mat;
}


}
#endif
