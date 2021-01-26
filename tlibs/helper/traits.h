/**
 * Custom type traits
 * @author Tobias Weber <tobias.weber@tum.de>
 * @date 19-nov-2014
 * @license GPLv2 or GPLv3
 */

#ifndef __MY_TRAITS_H__
#define __MY_TRAITS_H__

#include "boost_hacks.h"
#include <type_traits>


namespace tl {

// -----------------------------------------------------------------------------
template<class T, bool bScalar=std::is_scalar<T>::value>
struct underlying_value_type
{};

template<class T>
struct underlying_value_type<T, 1>
{
	using value_type = T;
};

template<class T>
struct underlying_value_type<T, 0>
{
	using value_type = typename underlying_value_type<
		typename T::value_type>::value_type;
};

template<class T>
using underlying_value_type_t =
	typename underlying_value_type<T, std::is_scalar<T>::value>::value_type;
// -----------------------------------------------------------------------------



// -----------------------------------------------------------------------------
}

#include <vector>
#include <array>
#include <list>
#include <initializer_list>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/math/quaternion.hpp>

namespace tl {

typedef std::integral_constant<int, 0> dim_0d_type;
typedef std::integral_constant<int, 1> dim_1d_type;
typedef std::integral_constant<int, 2> dim_2d_type;

template<class> struct get_type_dim : dim_0d_type {};

template<class... PARAMS> struct get_type_dim<std::vector<PARAMS...>> : dim_1d_type {};
template<class... PARAMS> struct get_type_dim<std::array<PARAMS...>> : dim_1d_type {};
template<class... PARAMS> struct get_type_dim<std::list<PARAMS...>> : dim_1d_type {};
template<class... PARAMS> struct get_type_dim<boost::numeric::ublas::vector<PARAMS...>> : dim_1d_type {};
template<class... PARAMS> struct get_type_dim<std::initializer_list<PARAMS...>> : dim_1d_type {};

template<class... PARAMS> struct get_type_dim<boost::numeric::ublas::matrix<PARAMS...>> : dim_2d_type {};



enum class LinalgType : short
{
	UNKNOWN,
	VECTOR,
	MATRIX,
	QUATERNION,
	REAL,
	COMPLEX
};

template<LinalgType val> struct linalg_type { static constexpr LinalgType value = val; };

template<class> struct get_linalg_type : linalg_type<LinalgType::UNKNOWN> {};
template<class... PARAMS> struct get_linalg_type<std::vector<PARAMS...>> : linalg_type<LinalgType::VECTOR> {};
template<class... PARAMS> struct get_linalg_type<boost::numeric::ublas::vector<PARAMS...>> : linalg_type<LinalgType::VECTOR> {};
template<class... PARAMS> struct get_linalg_type<boost::numeric::ublas::matrix<PARAMS...>> : linalg_type<LinalgType::MATRIX> {};
template<class... PARAMS> struct get_linalg_type<boost::math::quaternion<PARAMS...>> : linalg_type<LinalgType::QUATERNION> {};
template<class... PARAMS> struct get_linalg_type<std::complex<PARAMS...>> : linalg_type<LinalgType::COMPLEX> {};
template<> struct get_linalg_type<double> : linalg_type<LinalgType::REAL> {};
template<> struct get_linalg_type<float> : linalg_type<LinalgType::REAL> {};
// -----------------------------------------------------------------------------



// -----------------------------------------------------------------------------
}

#include <boost/units/dimensionless_quantity.hpp>

namespace tl {
enum class ScalarType : short
{
	TRIVIAL,
	DIMENSIONLESS,
	QUANTITY
};

template<ScalarType val> struct scalar_type { static constexpr ScalarType value = val; };

template<class T> struct get_scalar_type : scalar_type<ScalarType::TRIVIAL> { using value_type = T; };
template<class Sys, class T> struct get_scalar_type<boost::units::dimensionless_quantity<Sys, T>> : scalar_type<ScalarType::DIMENSIONLESS> { using value_type = T; };
template<class Sys, class T> struct get_scalar_type<boost::units::quantity<Sys, T>> : scalar_type<ScalarType::QUANTITY> 
{ using value_type = typename boost::units::quantity<Sys, T>::value_type; };
// -----------------------------------------------------------------------------




// -----------------------------------------------------------------------------
template<class T>
struct remove_constref
{
	typedef typename std::remove_const<
		typename std::remove_reference<T>::type
			>::type type;
};

// like C++14 style
template<class T>
using remove_constref_t = typename remove_constref<T>::type;
// -----------------------------------------------------------------------------



// -----------------------------------------------------------------------------
}

#include <utility>

namespace tl{

#if __cplusplus < 201402L
	template<class T, T... idx>
	struct integer_sequence
	{
		using value_type = T;
		static constexpr std::size_t size() { return sizeof...(idx); }
		static constexpr std::array<value_type, size()> value = {{idx...}};
	};

	template<class T, T NUM, typename _dummy = void> struct _make_integer_sequence;

	// TODO: make recursive
	template<class T, T NUM>
	struct _make_integer_sequence<T, NUM, typename std::enable_if<NUM==16>::type>
	{ using value_type = integer_sequence<T, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15>; };
	template<class T, T NUM>
	struct _make_integer_sequence<T, NUM, typename std::enable_if<NUM==15>::type>
	{ using value_type = integer_sequence<T, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14>; };
	template<class T, T NUM>
	struct _make_integer_sequence<T, NUM, typename std::enable_if<NUM==14>::type>
	{ using value_type = integer_sequence<T, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13>; };
	template<class T, T NUM>
	struct _make_integer_sequence<T, NUM, typename std::enable_if<NUM==13>::type>
	{ using value_type = integer_sequence<T, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12>; };
	template<class T, T NUM>
	struct _make_integer_sequence<T, NUM, typename std::enable_if<NUM==12>::type>
	{ using value_type = integer_sequence<T, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11>; };
	template<class T, T NUM>
	struct _make_integer_sequence<T, NUM, typename std::enable_if<NUM==11>::type>
	{ using value_type = integer_sequence<T, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10>; };
	template<class T, T NUM>
	struct _make_integer_sequence<T, NUM, typename std::enable_if<NUM==10>::type>
	{ using value_type = integer_sequence<T, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9>; };
	template<class T, T NUM>
	struct _make_integer_sequence<T, NUM, typename std::enable_if<NUM==9>::type>
	{ using value_type = integer_sequence<T, 0, 1, 2, 3, 4, 5, 6, 7, 8>; };
	template<class T, T NUM>
	struct _make_integer_sequence<T, NUM, typename std::enable_if<NUM==8>::type>
	{ using value_type = integer_sequence<T, 0, 1, 2, 3, 4, 5, 6, 7>; };
	template<class T, T NUM>
	struct _make_integer_sequence<T, NUM, typename std::enable_if<NUM==7>::type>
	{ using value_type = integer_sequence<T, 0, 1, 2, 3, 4, 5, 6>; };
	template<class T, T NUM>
	struct _make_integer_sequence<T, NUM, typename std::enable_if<NUM==6>::type>
	{ using value_type = integer_sequence<T, 0, 1, 2, 3, 4, 5>; };
	template<class T, T NUM>
	struct _make_integer_sequence<T, NUM, typename std::enable_if<NUM==5>::type>
	{ using value_type = integer_sequence<T, 0, 1, 2, 3, 4>; };
	template<class T, T NUM>
	struct _make_integer_sequence<T, NUM, typename std::enable_if<NUM==4>::type>
	{ using value_type = integer_sequence<T, 0, 1, 2, 3>; };
	template<class T, T NUM>
	struct _make_integer_sequence<T, NUM, typename std::enable_if<NUM==3>::type>
	{ using value_type = integer_sequence<T, 0, 1, 2>; };
	template<class T, T NUM>
	struct _make_integer_sequence<T, NUM, typename std::enable_if<NUM==2>::type>
	{ using value_type = integer_sequence<T, 0, 1>; };
	template<class T, T NUM>
	struct _make_integer_sequence<T, NUM, typename std::enable_if<NUM==1>::type>
	{ using value_type = integer_sequence<T, 0>; };
	template<class T, T NUM>
	struct _make_integer_sequence<T, NUM, typename std::enable_if<NUM==0>::type>
	{ using value_type = integer_sequence<T>; };

	template<class T, T NUM>
	using make_integer_sequence = typename _make_integer_sequence<T, NUM>::value_type;

#else
	// simply use C++14 standard implementation

	template<class T, T... idx>
	using integer_sequence = std::integer_sequence<T, idx...>;
	template<class T, T NUM>
	using make_integer_sequence = std::make_integer_sequence<T, NUM>;

#endif


/**
 * function call implementation
 */
template<class t_func,
	class t_arg = double, template<class ...> class t_cont = std::vector,
	std::size_t... idx>
t_arg _call_impl(t_func func, const t_cont<t_arg>& args,
	const /*std::*/integer_sequence<std::size_t, idx...>&)
{
	return func(args[idx]...);
}

/**
 * function call implementation (specialisation for std::array)
 */
template<class t_func, class t_arg = double, std::size_t... idx>
t_arg _call_impl(t_func func, const std::array<t_arg, sizeof...(idx)>& args,
	const /*std::*/integer_sequence<std::size_t, idx...>&)
{
	return func(args[idx]...);
}


/**
 * call a function with the args from an STL container
 */
template<std::size_t iNumArgs, class t_func,
	class t_arg = double, template<class ...> class t_cont = std::vector>
t_arg call(t_func func, const t_cont<t_arg>& args)
{
	using t_seq = /*std::*/make_integer_sequence<std::size_t, iNumArgs>;
	return _call_impl<t_func, t_arg, t_cont>(func, args, t_seq());
}

/**
 * call a function with the args from a std::array
 */
template<std::size_t iNumArgs, class t_func, class t_arg = double>
t_arg call(t_func func, const std::array<t_arg, iNumArgs>& args)
{
	using t_seq = /*std::*/make_integer_sequence<std::size_t, iNumArgs>;
	return _call_impl<t_func, t_arg>(func, args, t_seq());
}


// -----------------------------------------------------------------------------


template<typename t_arg, std::size_t ...idx>
using _t_fkt_vararg_impl = t_arg(*)(
	typename std::remove_reference<
		decltype(std::declval<t_arg*>()[idx])
	>::type...);

template<typename t_arg, std::size_t ...idx>
static _t_fkt_vararg_impl<t_arg, idx...>
_tstfkt_vararg(const tl::integer_sequence<std::size_t, idx...>&)
{ return nullptr; /* not interested in return value, only its type */ }


/**
 * constructs a function type with 'iNumArgs' arguments: t_arg (*) (t_arg, t_arg, ...)
 */
template<typename t_arg, std::size_t iNumArgs>
using t_fkt_vararg = decltype(
	_tstfkt_vararg<t_arg>(
		make_integer_sequence<std::size_t, iNumArgs>()));


// -----------------------------------------------------------------------------

}
#endif
