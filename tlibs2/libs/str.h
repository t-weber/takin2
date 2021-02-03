/**
 * tlibs2
 * string library
 * @author Tobias Weber <tobias.weber@tum.de>, <tweber@ill.fr>
 * @date 2013-2020
 * @license GPLv3, see 'LICENSE' file
 * @desc Forked on 7-Nov-2018 from my privately and TUM-PhD-developed "tlibs" project (https://github.com/t-weber/tlibs).
 */

#ifndef __TLIBS2_STRINGS__
#define __TLIBS2_STRINGS__

#include <string>
#include <cstring>
#include <iostream>
#include <sstream>
#include <locale>
#include <limits>
#include <map>
#include <utility>
#include <algorithm>
#include <type_traits>
#include <cctype>
#include <cwctype>
#include <unordered_map>

#include <boost/tokenizer.hpp>
#include <boost/algorithm/string.hpp>

#include "log.h"
#include "expr.h"


namespace algo = boost::algorithm;


namespace tl2 {

template<class t_str/*=std::string*/, class t_val/*=double*/>
std::pair<bool, t_val> eval_expr(const t_str& str) noexcept;


// -----------------------------------------------------------------------------



static inline std::wstring str_to_wstr(const std::string& str)
{
	return std::wstring(str.begin(), str.end());
}


static inline std::string wstr_to_str(const std::wstring& str)
{
	return std::string(str.begin(), str.end());
}


// overloaded in case the string is already of correct type
static inline const std::wstring& str_to_wstr(const std::wstring& str) { return str; }
static inline const std::string& wstr_to_str(const std::string& str) { return str; }


// -----------------------------------------------------------------------------


template<class t_str=std::string>
t_str str_to_upper(const t_str& str)
{
	t_str strOut;
	std::locale loc;

	strOut = boost::to_upper_copy(str, loc);
	return strOut;
}


template<class t_str=std::string>
t_str str_to_lower(const t_str& str)
{
	t_str strLower;
	std::locale loc;

	strLower = boost::to_lower_copy(str, loc);
	return strLower;
}


// -----------------------------------------------------------------------------


template<class t_str=std::string>
t_str get_file_noext(const t_str& str, bool bToLower=0)
{
	std::size_t iPos = str.find_last_of('.');

	if(iPos == t_str::npos)
		return str;

	t_str strRet = str.substr(0, iPos);
	if(bToLower)
		strRet = str_to_lower(strRet);

	return strRet;
}


template<class t_str=std::string>
t_str get_fileext(const t_str& str, bool bToLower=0)
{
	std::size_t iPos = str.find_last_of('.');

	if(iPos == t_str::npos)
		return t_str();

	t_str strRet = str.substr(iPos+1);
	if(bToLower)
		strRet = str_to_lower(strRet);

	return strRet;
}


/**
 *  e.g. returns "tof" for "123.tof.bz2"
 */
template<class t_str=std::string>
t_str get_fileext2(const t_str& str, bool bToLower=0)
{
	std::size_t iPos = str.find_last_of('.');
	if(iPos == t_str::npos || iPos == 0)
		return t_str();

	t_str strFile = str.substr(0, iPos);
	return get_fileext(strFile, bToLower);
}


/**
 * e.g. returns "tof" for "123.tof.bz2" and for "123.tof"
 */
template<class t_str=std::string>
t_str get_fileext_nocomp(const t_str& str, bool bToLower=0)
{
	std::size_t iCnt = std::count(str.begin(), str.end(), '.');
	if(iCnt==0)
		return t_str();
	else if(iCnt==1)
		return get_fileext(str, bToLower);
	else
		return get_fileext2(str, bToLower);
}


template<class t_str=std::string>
t_str get_dir(const t_str& str, bool bToLower=0)
{
	using t_ch = typename t_str::value_type;
	std::size_t iPos = std::string::npos;

	if constexpr(std::is_same_v<t_ch, char>)
		str.find_last_of("\\/");
	else if constexpr(std::is_same_v<t_ch, wchar_t>)
		str.find_last_of(L"\\/");

	if(iPos == t_str::npos)
		return t_str();

	t_str strRet = str.substr(0, iPos);
	if(bToLower)
		strRet = str_to_lower(strRet);

	return strRet;
}


template<class t_str=std::string>
t_str get_file_nodir(const t_str& str, bool bToLower=0)
{
	using t_ch = typename t_str::value_type;
	std::size_t iPos = std::string::npos;

	if constexpr(std::is_same_v<t_ch, char>)
		str.find_last_of("\\/");
	else if constexpr(std::is_same_v<t_ch, wchar_t>)
		str.find_last_of(L"\\/");

	if(iPos == t_str::npos)
		return t_str();

	t_str strRet = str.substr(iPos+1);
	if(bToLower)
		strRet = str_to_lower(strRet);

	return strRet;
}


// -----------------------------------------------------------------------------



template<class t_str=std::string>
bool str_is_equal(const t_str& str0, const t_str& str1, bool bCase=0)
{
	if(bCase)
		return algo::equals(str0, str1, algo::is_equal());
	else
		return algo::equals(str0, str1, algo::is_iequal(std::locale()));
}


template<class t_str=std::string>
bool str_is_equal_to_either(const t_str& str0,
	const std::initializer_list<t_str>& lststr1, bool bCase=0)
{
	for(const t_str& str1 : lststr1)
		if(str_is_equal<t_str>(str0, str1, bCase))
			return true;
	return false;
}


template<class t_str=std::string>
bool str_contains(const t_str& str, const t_str& strSub, bool bCase=0)
{
	if(bCase)
		return algo::contains(str, strSub, algo::is_equal());
	else
		return algo::contains(str, strSub, algo::is_iequal(std::locale()));
}


// -----------------------------------------------------------------------------


template<class t_str=std::string>
void trim(t_str& str)
{
	using t_char = typename t_str::value_type;

	boost::trim_if(str, [](t_char c) -> bool
	{
		if constexpr(std::is_same_v<t_char, char>)
			return t_str(" \t\r").find(c) != t_str::npos;
		else if constexpr(std::is_same_v<t_char, wchar_t>)
			return t_str(L" \t\r").find(c) != t_str::npos;

		return false;
	});
}


template<class t_str=std::string>
t_str trimmed(const t_str& str)
{
	t_str strret = str;
	trim(strret);
	return strret;
}


// -----------------------------------------------------------------------------



/**
 * removes all occurrences of a char in a string
 */
template<class t_str=std::string>
t_str remove_char(const t_str& str, typename t_str::value_type ch)
{
	t_str strRet;

	for(typename t_str::value_type c : str)
		if(c != ch)
			strRet.push_back(c);

	return strRet;
}


/**
 * removes all occurrences of specified chars in a string
 */
template<class t_str=std::string>
t_str remove_chars(const t_str& str, const t_str& chs)
{
	t_str strRet;

	for(typename t_str::value_type c : str)
	{
		bool bRemove = 0;
		for(typename t_str::value_type ch : chs)
		{
			if(c == ch)
			{
				bRemove = 1;
				break;
			}
		}

		if(!bRemove)
			strRet.push_back(c);
	}

	return strRet;
}


/**
 * Removes substring between strStart and strEnd
 * @return Number of removed substrings
 */
template<class t_str = std::string>
std::size_t string_rm(t_str& str, const t_str& strStart, const t_str& strEnd)
{
	std::size_t iNumFound = 0;

	while(1)
	{
		std::size_t iStart = str.find(strStart);
		std::size_t iEnd = str.find(strEnd);

		if(iStart == t_str::npos || iEnd == t_str::npos)
			break;
		if(iStart >= iEnd)
			break;

		str.erase(iStart, iEnd-iStart+strEnd.length());
		++iNumFound;
	}

	return iNumFound;
}


// -----------------------------------------------------------------------------


template<class t_str=std::string>
t_str insert_before(const t_str& str, const t_str& strChar, const t_str& strInsert)
{
	std::size_t pos = str.find(strChar);
	if(pos==t_str::npos)
		return str;

	t_str strRet = str;
	strRet.insert(pos, strInsert);

	return strRet;
}


template<class t_str=std::string>
bool begins_with(const t_str& str, const t_str& strBeg, bool bCase=1)
{
	if(bCase)
		return algo::starts_with(str, strBeg, algo::is_equal());
	else
		return algo::starts_with(str, strBeg, algo::is_iequal(std::locale()));
}


template<class t_str=std::string>
bool ends_with(const t_str& str, const t_str& strEnd, bool bCase=1)
{
	if(bCase)
		return algo::ends_with(str, strEnd, algo::is_equal());
	else
		return algo::ends_with(str, strEnd, algo::is_iequal(std::locale()));
}


// -----------------------------------------------------------------------------


template<class t_str=std::string>
std::pair<t_str, t_str>
split_first(const t_str& str, const t_str& strSep, bool bTrim=0, bool bSeq=0)
{
	t_str str1, str2;

	std::size_t iLenTok = bSeq ? strSep.length() : 1;
	std::size_t ipos = bSeq ? str.find(strSep) : str.find_first_of(strSep);

	if(ipos != t_str::npos)
	{
		str1 = str.substr(0, ipos);
		if(ipos+iLenTok < str.length())
			str2 = str.substr(ipos+iLenTok, t_str::npos);
	}

	if(bTrim)
	{
		trim(str1);
		trim(str2);
	}

	return std::make_pair(str1, str2);
}


/**
 * get string between strSep1 and strSep2
 */
template<class t_str=std::string>
t_str str_between(const t_str& str, const t_str& strSep1, const t_str& strSep2,
	bool bTrim=0, bool bSeq=0)
{
	t_str str1, str2;
	std::tie(str1, str2) = split_first<t_str>(str, strSep1, bTrim, bSeq);
	if(str2 == "") return t_str("");

	std::tie(str1, str2) = split_first<t_str>(str2, strSep2, bTrim, bSeq);
	return str1;
}


// ----------------------------------------------------------------------------


template<typename T, class t_str=std::string, bool bTIsStr=0>
struct _str_to_var_impl;


template<typename T, class t_str>
struct _str_to_var_impl<T, t_str, 1>
{
	inline const T& operator()(const t_str& str) const
	{
		return str;
	}
};


template<typename T, class t_str>
struct _str_to_var_impl<T, t_str, 0>
{
	inline T operator()(const t_str& str) const
	{
		if(!trimmed(str).length())
			return T();

		typedef typename t_str::value_type t_char;
		std::basic_istringstream<t_char> istr(str);

		T t;
		istr >> t;
		return t;
	}
};



/**
 * tokenises string on any of the chars in strDelim
 */
template<class T, class t_str=std::string, class t_cont=std::vector<T>>
void get_tokens(const t_str& str, const t_str& strDelim, t_cont& vecRet)
{
	using t_char = typename t_str::value_type;
	using t_tokeniser = boost::tokenizer<boost::char_separator<t_char>,
		typename t_str::const_iterator, t_str>;
	using t_tokiter = typename t_tokeniser::iterator;

	boost::char_separator<t_char> delim(strDelim.c_str());
	t_tokeniser tok(str, delim);

	for(t_tokiter iter=tok.begin(); iter!=tok.end(); ++iter)
	{
		vecRet.push_back(
			_str_to_var_impl<T, t_str,
			std::is_convertible<T, t_str>::value>()(*iter));
	}
}


/**
 * Tokenises string on strDelim
 */
template<class T, class t_str=std::string, template<class...> class t_cont=std::vector>
void get_tokens_seq(const t_str& str, const t_str& strDelim, t_cont<T>& vecRet, bool bCase=1)
{
	using t_char = typename t_str::value_type;

	std::locale loc;
	t_cont<t_str> vecStr;
	algo::iter_split(vecStr, str, algo::first_finder(strDelim,
		[bCase, &loc](t_char c1, t_char c2) -> bool
		{
			if(!bCase)
			{
				c1 = std::tolower(c1, loc);
				c2 = std::tolower(c2, loc);
			}

			return c1 == c2;
		}));

	for(const t_str& strTok : vecStr)
	{
		vecRet.push_back(
			_str_to_var_impl<T, t_str,
			std::is_convertible<T, t_str>::value>()(strTok));
	}
}


template<class T, class t_str=std::string, class t_cont=std::vector<T>>
bool parse_tokens(const t_str& str, const t_str& strDelim, t_cont& vecRet)
{
	std::vector<t_str> vecStrs;
	get_tokens<t_str, t_str, std::vector<t_str>>(str, strDelim, vecStrs);

	bool bOk = 1;
	for(const t_str& str : vecStrs)
	{
		std::pair<bool, T> pairResult = eval_expr<t_str, T>(str);
		vecRet.push_back(pairResult.second);
		if(!pairResult.first) bOk = 0;
	}

	return bOk;
}


template<typename T, class t_str=std::string>
T str_to_var_parse(const t_str& str)
{
	std::pair<bool, T> pairResult = eval_expr<t_str, T>(str);
	if(!pairResult.first)
		return T(0);
	return pairResult.second;
}


template<typename T, class t_str=std::string>
T str_to_var(const t_str& str)
{
	return _str_to_var_impl<T, t_str, std::is_convertible<T, t_str>::value>()(str);
}


// ----------------------------------------------------------------------------


template<class T, class t_ch,
	bool is_number_type=std::is_fundamental<T>::value>
struct _var_to_str_print_impl {};


template<class T, class t_ch> struct _var_to_str_print_impl<T, t_ch, false>
{
	void operator()(std::basic_ostream<t_ch>& ostr, const T& t) { ostr << t; }
};


template<class T, class t_ch> struct _var_to_str_print_impl<T, t_ch, true>
{
	void operator()(std::basic_ostream<t_ch>& ostr, const T& t)
	{
		// prevents printing "-0"
		T t0 = t;
		if(t0==T(-0)) t0 = T(0);

		ostr << t0;
	}
};


template<typename T, class t_str=std::string>
struct _var_to_str_impl
{
	t_str operator()(const T& t,
		std::streamsize iPrec = std::numeric_limits<T>::max_digits10,
		int iGroup=-1)
	{
		//if(std::is_convertible<T, t_str>::value)
		//	return *reinterpret_cast<const t_str*>(&t);

		typedef typename t_str::value_type t_char;

		std::basic_ostringstream<t_char> ostr;
		ostr.precision(iPrec);


		class Sep : public std::numpunct<t_char>
		{
		public:
			Sep() : std::numpunct<t_char>(1) {}
			~Sep() { /*std::cout << "~Sep();" << std::endl;*/ }
		protected:
			virtual t_char do_thousands_sep() const override { return ' ';}
			virtual std::string do_grouping() const override { return "\3"; }
		};
		Sep *pSep = nullptr;


		if(iGroup > 0)
		{
			pSep = new Sep();
			ostr.imbue(std::locale(ostr.getloc(), pSep));
		}

		_var_to_str_print_impl<T, t_char> pr;
		pr(ostr, t);
		t_str str = ostr.str();

		if(pSep)
		{
			ostr.imbue(std::locale());
			delete pSep;
		}
		return str;
	}
};


template<class t_str>
struct _var_to_str_impl<t_str, t_str>
{
	const t_str& operator()(const t_str& tstr, std::streamsize /*iPrec=10*/, int /*iGroup=-1*/)
	{
		return tstr;
	}

	t_str operator()(const typename t_str::value_type* pc, std::streamsize /*iPrec=10*/, int /*iGroup=-1*/)
	{
		return t_str(pc);
	}
};


template<typename T, class t_str=std::string>
t_str var_to_str(const T& t,
	std::streamsize iPrec = std::numeric_limits<T>::max_digits10-1,
	int iGroup = -1)
{
	_var_to_str_impl<T, t_str> _impl;
	return _impl(t, iPrec, iGroup);
}


/**
 * converts a container (e.g. a vector) to a string
 */
template<class t_cont, class t_str=std::string>
t_str cont_to_str(const t_cont& cont, const char* pcDelim=",",
	std::streamsize iPrec = std::numeric_limits<typename t_cont::value_type>::max_digits10-1)
{
	using t_elem = typename t_cont::value_type;

	t_str str;

	for(typename t_cont::const_iterator iter=cont.begin(); iter!=cont.end(); ++iter)
	{
		const t_elem& elem = *iter;

		str += var_to_str<t_elem, t_str>(elem, iPrec);
		if(iter+1 != cont.end())
			str += pcDelim;
	}
	return str;
}

// ----------------------------------------------------------------------------


template<typename t_char=char>
bool skip_after_line(std::basic_istream<t_char>& istr,
	const std::basic_string<t_char>& strLineBegin,
	bool bTrim=true, bool bCase=0)
{
	while(!istr.eof())
	{
		std::basic_string<t_char> strLine;
		std::getline(istr, strLine);
		if(bTrim)
			trim(strLine);

		if(strLine.size() < strLineBegin.size())
			continue;

		std::basic_string<t_char> strSub = strLine.substr(0, strLineBegin.size());

		if(str_is_equal<std::basic_string<t_char>>(strSub, strLineBegin, bCase))
			return true;
	}
	return false;
}


template<typename t_char=char>
void skip_after_char(std::basic_istream<t_char>& istr, t_char ch, bool bCase=0)
{
	std::locale loc;
	if(!bCase) ch = std::tolower(ch, loc);

	while(!istr.eof())
	{
		t_char c;
		istr.get(c);

		if(!bCase) c = std::tolower(c, loc);

		if(c == ch)
			break;
	}
}


// ----------------------------------------------------------------------------


/**
 * functions working on chars
 */
template<class t_ch, typename=void> struct char_funcs {};


/**
 * specialisation for char
 */
template<class t_ch>
struct char_funcs<t_ch, typename std::enable_if<std::is_same<t_ch, char>::value>::type>
{
	static bool is_digit(t_ch c) { return std::isdigit(c); }
};


/**
 * specialisation for wchar_t
 */
template<class t_ch>
struct char_funcs<t_ch, typename std::enable_if<std::is_same<t_ch, wchar_t>::value>::type>
{
	static bool is_digit(t_ch c) { return std::iswdigit(c); }
};


// ----------------------------------------------------------------------------


/**
 * tests if a string consists entirely of numbers
 */
template<class t_str = std::string>
bool str_is_digits(const t_str& str)
{
	using t_ch = typename t_str::value_type;
	using t_fkt = char_funcs<t_ch>;

	bool bAllNums = std::all_of(str.begin(), str.end(),
		[](t_ch c) -> bool
		{
			return t_fkt::is_digit(c);
		});
	return bAllNums;
}



// ----------------------------------------------------------------------------



template<class t_str=std::string, class t_cont=std::vector<double>>
t_cont get_py_array(const t_str& str)
{
	typedef typename t_cont::value_type t_elems;
	t_cont vecArr;

	std::size_t iStart = str.find('[');
	std::size_t iEnd = str.find(']');

	// search for list instead
	if(iStart==t_str::npos || iEnd==t_str::npos)
	{
		iStart = str.find('(');
		iEnd = str.find(')');
	}

	// invalid array
	if(iStart==t_str::npos || iEnd==t_str::npos || iEnd<iStart)
		return vecArr;

	t_str strArr = str.substr(iStart+1, iEnd-iStart-1);
	get_tokens<t_elems, t_str>(strArr, ",", vecArr);

	return vecArr;
}


template<class t_str=std::string>
t_str get_py_string(const t_str& str)
{
	std::size_t iStart = str.find_first_of("\'\"");
	std::size_t iEnd = str.find_last_of("\'\"");

	// invalid string
	if(iStart==t_str::npos || iEnd==t_str::npos || iEnd<iStart)
		return "";

	return str.substr(iStart+1, iEnd-iStart-1);
}



// ----------------------------------------------------------------------------



template<class t_str/*=std::string*/, class t_val/*=double*/>
std::pair<bool, t_val> eval_expr(const t_str& str) noexcept
{
	if(trimmed(str).length() == 0)
		return std::make_pair(true, t_val(0));

	try
	{
		ExprParser<t_val> parser;
		bool ok = parser.parse(wstr_to_str(str));
		t_val valRes = parser.eval();
		return std::make_pair(ok, valRes);
	}
	catch(const std::exception& ex)
	{
		log_err("Parsing failed with error: ", ex.what(), ".");
		return std::make_pair(false, t_val(0));
	}
}

}
#endif
