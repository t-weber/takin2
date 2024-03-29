/**
 * Wrapper for clipper spacegroups
 * @author Tobias Weber <tobias.weber@tum.de>
 * @date oct-2015
 * @license GPLv2
 *
 * ----------------------------------------------------------------------------
 * Takin (inelastic neutron scattering software package)
 * Copyright (C) 2017-2021  Tobias WEBER (Institut Laue-Langevin (ILL),
 *                          Grenoble, France).
 * Copyright (C) 2013-2017  Tobias WEBER (Technische Universitaet Muenchen
 *                          (TUM), Garching, Germany).
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; version 2 of the License.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along
 * with this program; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 * ----------------------------------------------------------------------------
 */

#ifndef __XTL_SGCLP_H__
#define __XTL_SGCLP_H__

#include <map>
#include <clipper/clipper.h>

#include "crystalsys.h"
#include "tlibs/string/string.h"
#include "sghelper.h"


namespace xtl {

namespace ublas = boost::numeric::ublas;


class SpaceGroupClp
{
protected:
	clipper::Spacegroup *m_psg = nullptr;
	std::string m_strName;
	std::string m_strLaue;
	CrystalSystem m_crystalsys;

public:
	SpaceGroupClp(unsigned int iNum);
	SpaceGroupClp(const SpaceGroupClp& sg);
	SpaceGroupClp(SpaceGroupClp&& sg);
	SpaceGroupClp() = delete;
	~SpaceGroupClp();

	bool HasReflection(int h, int k, int l) const;
	bool HasReflection2(int h, int k, int l) const;

	void SetName(const std::string& str) { m_strName = str; }
	const std::string& GetName() const { return m_strName; }

	const std::string& GetLaueGroup() const { return m_strLaue; }

	CrystalSystem GetCrystalSystem() const { return m_crystalsys; }
	const char* GetCrystalSystemName() const { return get_crystal_system_name(m_crystalsys); }

	void GetSymTrafos(std::vector<ublas::matrix<double>>& vecTrafos) const;
	void GetInvertingSymTrafos(std::vector<ublas::matrix<double>>& vecTrafos) const;
	void GetPrimitiveSymTrafos(std::vector<ublas::matrix<double>>& vecTrafos) const;
	void GetCenteringSymTrafos(std::vector<ublas::matrix<double>>& vecTrafos) const;
};


typedef std::map<std::string, SpaceGroupClp> t_mapSpaceGroups;
extern const t_mapSpaceGroups* get_space_groups();
extern void init_space_groups();


template<class T=double, class t_clp=clipper::ftype>
ublas::matrix<T> symop_to_matrix(const clipper::Symop& symop)
{
	const clipper::Mat33<t_clp>& mat = symop.rot();
	const clipper::Vec3<t_clp>& vec = symop.trn();

	ublas::matrix<T> matRet(4,4);
	for(int i=0; i<3; ++i)
		for(int j=0; j<3; ++j)
			matRet(i,j) = T(mat(i,j));

	for(int i=0; i<3; ++i)
	{
		matRet(i,3) = T(vec[i]);
		matRet(3,i) = 0.;
	}

	matRet(3,3) = 1.;
	return matRet;
}


enum class SymTrafoType
{
	ALL,

	PRIMITIVE,
	INVERTING,
	CENTERING
};

template<class T=clipper::ftype>
void get_symtrafos(const clipper::Spacegroup& sg,
	std::vector<ublas::matrix<T>>& vecTrafos, SymTrafoType ty=SymTrafoType::ALL)
{
	int iNumSymOps = 0;

	if(ty == SymTrafoType::ALL)
		iNumSymOps = sg.num_symops();
	else if(ty == SymTrafoType::INVERTING)
		iNumSymOps = sg.num_inversion_symops();
	else if(ty == SymTrafoType::PRIMITIVE)
		iNumSymOps = sg.num_primitive_symops();
	else if(ty == SymTrafoType::CENTERING)
		iNumSymOps = sg.num_centering_symops();

	vecTrafos.clear();
	vecTrafos.reserve(iNumSymOps);

	for(int iSymOp=0; iSymOp<iNumSymOps; ++iSymOp)
	{
		const clipper::Symop* symop = nullptr;

		if(ty == SymTrafoType::ALL)
			symop = &sg.symop(iSymOp);
		else if(ty == SymTrafoType::INVERTING)
			symop = &sg.inversion_symop(iSymOp);
		else if(ty == SymTrafoType::PRIMITIVE)
			symop = &sg.primitive_symop(iSymOp);
		else if(ty == SymTrafoType::CENTERING)
			symop = &sg.centering_symop(iSymOp);

		ublas::matrix<T> mat = symop_to_matrix<T>(*symop);
		vecTrafos.push_back(mat);
	}
}

}
#endif
