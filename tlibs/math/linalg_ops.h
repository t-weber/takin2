/**
 * Linalg operators
 * @author Tobias Weber <tobias.weber@tum.de>
 * @date apr-2015
 * @license GPLv2 or GPLv3
 */

#ifndef __TL_LINALG_OPS_H__
#define __TL_LINALG_OPS_H__

#include <boost/numeric/ublas/vector_expression.hpp>
#include <boost/numeric/ublas/matrix_expression.hpp>

#include "math.h"
#include "linalg.h"
#include "../helper/exception.h"
#include "../helper/traits.h"

namespace tl {

namespace ublas = boost::numeric::ublas;

template<class T1, class T2,
	LinalgType ty1=get_linalg_type<T1>::value,
	LinalgType ty2=get_linalg_type<T2>::value>
struct linalg_mult_op_impl;


// vec * vec
template<class T1, class T2>
struct linalg_mult_op_impl<T1, T2,
	LinalgType::VECTOR, LinalgType::VECTOR>
{
	typedef typename T1::value_type ret_type;

	ret_type operator()(const T1& vec1, const T2& vec2) const
	{
		//typedef typename T1::value_type REAL;
		return inner(vec1, vec2);
	}
};


// mat * mat
template<class T1, class T2>
struct linalg_mult_op_impl<T1, T2,
	LinalgType::MATRIX, LinalgType::MATRIX>
{
	typedef T1 ret_type;

	ret_type operator()(const T1& mat1, const T2& mat2) const
	{
		return prod_mm(mat1, mat2);
	}
};


// mat * vec
template<class T1, class T2>
struct linalg_mult_op_impl<T1, T2,
	LinalgType::MATRIX, LinalgType::VECTOR>
{
	typedef T2 ret_type;

	ret_type operator()(const T1& mat, const T2& vec) const
	{
		return prod_mv(mat, vec);
	}
};



template<class T1, class T2>
typename tl::linalg_mult_op_impl<T1, T2>::ret_type mult(
	const typename std::enable_if<tl::get_linalg_type<T1>::value != tl::LinalgType::UNKNOWN, T1>::type& t1,
	const typename std::enable_if<tl::get_linalg_type<T1>::value != tl::LinalgType::UNKNOWN, T2>::type& t2)
{
	return tl::linalg_mult_op_impl<T1, T2>()(t1, t2);
}


// ----------------------------------------------------------------------------


template<class T1, class T2, class EPS,
	LinalgType ty1 = get_linalg_type<T1>::value,
	LinalgType ty2 = get_linalg_type<T2>::value>
struct linalg_equ_op_impl;


template<class T1, class T2, class EPS>
struct linalg_equ_op_impl<T1, T2, EPS,
	LinalgType::MATRIX, LinalgType::MATRIX>
{
	bool operator()(const T1& mat1, const T2& mat2, EPS eps) const
	{
		return mat_equal(mat1, mat2, eps);
	}
};


template<class T1, class T2, class EPS>
struct linalg_equ_op_impl<T1, T2, EPS,
	LinalgType::VECTOR, LinalgType::VECTOR>
{
	bool operator()(const T1& vec1, const T2& vec2, EPS eps) const
	{
		return vec_equal(vec1, vec2, eps);
	}
};


template<class T1, class T2>
bool equ(
	const typename std::enable_if<tl::get_linalg_type<T1>::value != tl::LinalgType::UNKNOWN, T1>::type& t1,
	const typename std::enable_if<tl::get_linalg_type<T2>::value != tl::LinalgType::UNKNOWN, T2>::type& t2,
	typename tl::_get_epsilon_impl<T1>::t_eps eps = tl::get_epsilon<T1>())
{
	return tl::linalg_equ_op_impl<T1, T2, decltype(eps)>()(t1, t2, eps);
}

}


#ifndef TLIBS_USE_GLOBAL_OPS
namespace tl_ops {
#endif

template<class T1, class T2>
typename tl::linalg_mult_op_impl<T1, T2>::ret_type operator*
	(const T1& t1, const T2& t2)
{
	return tl::mult<T1, T2>(t1, t2);
}


template<class T1, class T2>
bool operator==(const T1& t1, const T2& t2)
{
	typename tl::_get_epsilon_impl<T1>::t_eps eps = tl::get_epsilon<T1>();
	return tl::equ<T1, T2>(t1, t2, eps);
}

#ifndef TLIBS_USE_GLOBAL_OPS
}	// tl_ops
#endif

#endif
