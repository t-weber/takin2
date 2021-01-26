/**
 * property tree wrapper
 * @author Tobias Weber <tobias.weber@tum.de>
 * @date dec-2015
 * @license GPLv2 or GPLv3
 */

#ifndef __PROP_FILES_H__
#define __PROP_FILES_H__

#include <string>
#include <iostream>
#include <sstream>
#include <fstream>
#include <type_traits>
#include <map>
#include <algorithm>

#include <boost/version.hpp>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/xml_parser.hpp>
#include <boost/property_tree/json_parser.hpp>
#include <boost/property_tree/info_parser.hpp>
#include <boost/property_tree/ini_parser.hpp>
#include <boost/optional.hpp>

#include "../string/string.h"
#include "../helper/traits.h"
#include "../helper/misc.h"

#if !defined NO_STR_PARSER
	#include "../string/eval.h"
#endif

#if !defined NO_IOSTR
	#include "../file/comp.h"
#else
	enum class Compressor { INVALID };
#endif


namespace tl {

namespace prop = ::boost::property_tree;

enum class PropType { XML, JSON, INFO, INI, };

// ----------------------------------------------------------------------------
// string compare predicate
template<class t_str, bool bCaseSensitive> struct StringComparer {};

// honour case
template<class t_str> struct StringComparer<t_str, 1>
{
	bool operator()(const t_str& str1, const t_str& str2) const
	{
		return std::lexicographical_compare(str1.begin(), str1.end(),
			str2.begin(), str2.end());
	}
};

// ignore case
template<class t_str> struct StringComparer<t_str, 0>
{
	bool operator()(const t_str& _str1, const t_str& _str2) const
	{
		t_str str1 = str_to_lower<t_str>(_str1);
		t_str str2 = str_to_lower<t_str>(_str2);

		return std::lexicographical_compare(str1.begin(), str1.end(),
			str2.begin(), str2.end());
	}
};

// ----------------------------------------------------------------------------


template<class t_str = std::string>
boost::optional<prop::string_path<t_str, prop::id_translator<t_str>>>
get_prop_path(const std::string& _strAddr, const typename t_str::value_type& chSep)
{
	using t_ret = boost::optional<prop::string_path<t_str, prop::id_translator<t_str>>>;

	t_str strAddr = _strAddr;
	trim(strAddr);

	if(strAddr.length() == 0)
		return t_ret();

	if(strAddr[0] == chSep)
		strAddr = strAddr.substr(1);

	try
	{
		prop::string_path<t_str, prop::id_translator<t_str>> path(strAddr, chSep);
		return boost::optional<decltype(path)>(path);
	}
	catch(const prop::ptree_bad_path& ex) {}
	catch(const std::exception& ex) {}

	return t_ret();
}


// ----------------------------------------------------------------------------


template<class _t_str = std::string, bool bCaseSensitive=0>
class Prop
{
public:
	using t_str = _t_str;
	using t_ch = typename t_str::value_type;
	using t_prop = prop::basic_ptree<t_str, t_str, StringComparer<t_str, bCaseSensitive>>;
	using t_propval = typename t_prop::value_type;

protected:
	t_prop m_prop;
	t_ch m_chSep = '/';

public:
	Prop(t_ch chSep = '/') : m_chSep(chSep) {}
	Prop(const t_prop& prop, t_ch chSep='/') : m_prop(prop), m_chSep(chSep) {}
	Prop(t_prop&& prop, t_ch chSep='/') : m_prop(std::move(prop)), m_chSep(chSep) {}
	virtual ~Prop() = default;

	void SetSeparator(t_ch ch) { m_chSep = ch; }

	const t_prop& GetProp() const { return m_prop; }
	void SetProp(const t_prop& prop) { m_prop = prop; }

	bool Load(const t_ch* pcFile)
	{
		t_str strFile = pcFile;
		t_str strExt = str_to_lower<t_str>(get_fileext_nocomp<t_str>(strFile));

		if(strExt == "xml")
			return Load(pcFile, PropType::XML);
		else if(strExt == "json")
			return Load(pcFile, PropType::JSON);
		else if(strExt == "info")
			return Load(pcFile, PropType::INFO);
		else if(strExt == "ini")
			return Load(pcFile, PropType::INI);

		return false;
	}


	bool Load(const t_ch* pcFile, PropType ty)
	{
		std::basic_ifstream<t_ch> ifstr(pcFile);
		if(!ifstr) return false;

	#if !defined NO_IOSTR
		std::shared_ptr<std::basic_istream<t_ch>> ptrIstr
			= create_autodecomp_istream(ifstr);
		if(!ptrIstr)
			return false;
		std::basic_istream<t_ch>* pIfstr = ptrIstr.get();
	#else
		std::basic_istream<t_ch>* pIfstr = &ifstr;
	#endif

		if(!*pIfstr) return false;
		return Load(*pIfstr, ty);
	}


	bool Load(std::basic_istream<t_ch>& istr, PropType ty)
	{
		try
		{
			switch(ty)
			{
				case PropType::XML:
					prop::read_xml(istr, m_prop);
					break;
				case PropType::JSON:
					prop::read_json(istr, m_prop);
					break;
				case PropType::INFO:
					prop::read_info(istr, m_prop);
					break;
				case PropType::INI:
					prop::read_ini(istr, m_prop);
					break;
				default:
					return false;
			}
		}
		catch(const prop::file_parser_error& err)
		{
			return false;
		}
		return true;
	}


	bool Save(const t_ch* pcFile) const
	{
		t_str strFile = pcFile;
		t_str strExt = str_to_lower<t_str>(get_fileext_nocomp<t_str>(strFile));
		t_str strComp = str_to_lower<t_str>(get_fileext<t_str>(strFile));
		Compressor comp = Compressor::INVALID;
	#if !defined NO_IOSTR
		comp = comp_from_ext(strComp);
	#endif

		if(strExt == "xml")
			return Save(pcFile, PropType::XML, comp);
		else if(strExt == "json")
			return Save(pcFile, PropType::JSON, comp);
		else if(strExt == "info")
			return Save(pcFile, PropType::INFO, comp);
		else if(strExt == "ini")
			return Save(pcFile, PropType::INI, comp);

		return false;
	}


	bool Save(const t_ch* pcFile, PropType ty, Compressor comp=Compressor::INVALID) const
	{
		std::basic_ofstream<t_ch> ofstr(pcFile);
		if(!ofstr) return false;

	#if !defined NO_IOSTR
		std::shared_ptr<std::basic_ostream<t_ch>> ptrOstr
			= create_comp_ostream(ofstr, comp);
		if(!ptrOstr)
			return false;
		std::basic_ostream<t_ch>* pOfstr = ptrOstr.get();
	#else
		std::basic_ostream<t_ch>* pOfstr = &ofstr;
	#endif

		if(!*pOfstr) return false;
		return Save(*pOfstr, ty);
	}


	bool Save(std::basic_ostream<t_ch>& ofstr, PropType ty) const
	{
		#if BOOST_VERSION >= 105700
			using t_writer = t_str;
		#else
			using t_writer = t_ch;
		#endif

		try
		{
			switch(ty)
			{
				case PropType::XML:
					prop::write_xml(ofstr, m_prop, prop::xml_writer_settings<t_writer>('\t',1));
					break;
				case PropType::JSON:
					prop::write_json(ofstr, m_prop);
					break;
				case PropType::INFO:
					prop::write_info(ofstr, m_prop);
					break;
				case PropType::INI:
					prop::write_ini(ofstr, m_prop);
					break;
				default:
					return false;
			}
		}
		catch(const prop::file_parser_error& err)
		{
			return false;
		}

		return true;
	}


	bool Load(const t_str& strFile) { return Load(strFile.c_str()); }
	bool Load(const t_str& strFile, PropType ty) { return Load(strFile.c_str(), ty); }
	bool Save(const t_str& strFile) const { return Save(strFile.c_str()); }
	bool Save(const t_str& strFile, PropType ty, Compressor comp=Compressor::INVALID) const
	{ return Save(strFile.c_str(), ty, comp); }


	template<typename T>
	T Query(const t_str& strAddr, const T* pDef=nullptr, bool *pbOk=nullptr) const
	{
		T tOut;
		try
		{
			auto optPath = get_prop_path<t_str>(strAddr, m_chSep);
			if(!optPath) throw std::exception();

			tOut = tl::str_to_var<T, t_str>(m_prop.template get<t_str>(*optPath));
		}
		catch(const std::exception& ex)
		{
			if(pbOk) *pbOk = 0;
			if(pDef) return *pDef;
			return T();
		}

		// if T is a string type, trim it
		if(std::is_same<t_str, T>::value)
			trim(*reinterpret_cast<t_str*>(&tOut));

		if(pbOk) *pbOk = 1;
		return tOut;
	}


	template<typename T>
	T Query(const t_str& _strAddr, const T def, bool *pbOk=nullptr) const
	{
		return Query<T>(_strAddr, &def, pbOk);
	}


	template<typename T>
	boost::optional<T> QueryOpt(const t_str& strAddr) const
	{
		bool bOk = 0;
		T tVal = Query<T>(strAddr, nullptr, &bOk);
		return bOk ? boost::optional<T>(std::move(tVal)) : boost::optional<T>();
	}


#if !defined NO_STR_PARSER
	template<typename T>
	T QueryAndParse(const t_str& _strAddr, const T* pDef=nullptr, bool *pbOk=nullptr) const
	{
		bool bOk = 0;
		T def = pDef ? * pDef : T();
		t_str strExpr = Query<t_str>(_strAddr, nullptr, &bOk);

		if(pbOk) *pbOk = bOk;
		if(!bOk) return def;

		std::pair<bool, T> pairRes = eval_expr<t_str, T>(strExpr);
		if(!pairRes.first)
		{
			if(pbOk) *pbOk = 0;
			return def;
		}

		return pairRes.second;
	}

	template<typename T>
	T QueryAndParse(const t_str& _strAddr, const T def, bool *pbOk=nullptr) const
	{
		return QueryAndParse<T>(_strAddr, &def, pbOk);
	}

#else 	// simply call normal query function

	template<typename T>
	T Query(const t_str& _strAddr, const T* pDef=nullptr, bool *pbOk=nullptr) const
	{
		return Query<T>(_strAddr, pDef, pbOk);
	}

	template<typename T>
	T QueryAndParse(const t_str& _strAddr, const T def, bool *pbOk=nullptr) const
	{
		return Query<T>(_strAddr, def, pbOk);
	}

#endif


	/**
	 * get children to a node
	 */
	std::vector<t_propval> GetFullChildNodes(const t_str& strAddr) const
	{
		std::vector<t_propval> vecRet;
		try
		{
			auto optPath = get_prop_path<t_str>(strAddr, m_chSep);
			if(!optPath) throw std::exception();

			for(const auto &node : m_prop.get_child(*optPath))
				vecRet.push_back(node);
		}
		catch(const std::exception& ex) {}

		return vecRet;
	}


	/**
	 * get a list of children to a node
	 */
	std::vector<t_str> GetChildNodes(const t_str& strAddr) const
	{
		std::vector<t_str> vecRet;
		try
		{
			auto optPath = get_prop_path<t_str>(strAddr, m_chSep);
			if(!optPath) throw std::exception();

			for(const auto &node : m_prop.get_child(*optPath))
				vecRet.push_back(node.first);
		}
		catch(const std::exception& ex) {}

		return vecRet;
	}


	/**
	 * get a list of child values to a node
	 */
	template<class T = t_str>
	std::vector<T> GetChildValues(const t_str& strAddr) const
	{
		std::vector<T> vecRet;
		try
		{
			auto optPath = get_prop_path<t_str>(strAddr, m_chSep);
			if(!optPath) throw std::exception();

			for(const auto &node : m_prop.get_child(*optPath))
				vecRet.push_back(node.second.template get_value<T>());
		}
		catch(const std::exception& ex) {}

		return vecRet;
	}


	bool Exists(const t_str& strAddr) const
	{
		bool bOk = 0;
		t_str strQuery = Query<t_str>(strAddr, nullptr, &bOk);
		if(strQuery.length() == 0)
			bOk = 0;

		return bOk;
	}


	bool PathExists(const t_str& strAddr) const
	{
		return GetChildNodes(strAddr).size() != 0;
	}


	template<class T>
	void Add(const t_str& strKey, T&& tVal)
	{
		auto optPath = get_prop_path<t_str>(strKey, m_chSep);
		if(!optPath) return;

		if(std::is_convertible<T, t_str>::value)
		{
			m_prop.add(std::move(*optPath), std::forward<T>(tVal));
		}
		else
		{
			t_str strVal = var_to_str<T, t_str>(tVal);
			m_prop.add(std::move(*optPath), std::move(strVal));
		}
	}


	template<class t_map = std::map<t_str, t_str>>
	void Add(const t_map& map)
	{
		using t_key = typename t_map::key_type;
		using t_val = typename t_map::mapped_type;
		using t_pair = typename t_map::value_type;

		for(const t_pair& pair : map)
		{
			t_str strKey = var_to_str<t_key, t_str>(pair.first);
			t_str strVal = var_to_str<t_val, t_str>(pair.second);

			Add<t_str>(std::move(strKey), std::move(strVal));
		}
	}


	friend std::ostream& operator<<(std::ostream& ostr,
		const Prop<_t_str, bCaseSensitive>& prop)
	{
		return operator << (ostr, prop.m_prop);
	}
};


template<class t_str = std::string, bool bCaseSens = 0>
std::ostream& operator<<(std::ostream& ostr,
	const prop::basic_ptree<t_str, t_str, StringComparer<t_str, bCaseSens>>& prop)
{
	using t_prop = prop::basic_ptree<t_str, t_str, StringComparer<t_str, bCaseSens>>;

	for(typename t_prop::const_iterator iter = prop.begin();
		iter!= prop.end(); ++iter)
	{
		std::cout << iter->first << " " << iter->second.data() << " "
			<< iter->second << "\n";
	}
	return ostr;
}

}
#endif
