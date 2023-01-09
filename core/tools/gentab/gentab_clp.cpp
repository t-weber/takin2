/**
 * Creates needed tables
 * @author Tobias Weber <tobias.weber@tum.de>
 * @date nov-2015
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


// ----------------------------------------------------------------------------
// ugly, but can't be helped for the moment:
// directly link to the internal clipper coefficients table
// that lives in clipper/clipper/core/atomsf.cpp
namespace clipper { namespace data
{
	extern const struct SFData
	{
		const char atomname[8];
		const t_real a[5], c, b[5], d;  // d is always 0
	} sfdata[];

	const unsigned int numsfdata = 212;
}}
// ----------------------------------------------------------------------------

namespace dat = clipper::data;


// ============================================================================


bool gen_formfacts_clp()
{
	tl::Prop<std::string> prop;
	prop.SetSeparator('.');

	prop.Add("ffacts.source", "Form factor coefficients extracted from Clipper (by K. Cowtan).");
	prop.Add("ffacts.source_url", "http://www.ysbl.york.ac.uk/~cowtan/clipper/");
	prop.Add("ffacts.num_atoms", tl::var_to_str(dat::numsfdata));

	for(unsigned int iFF=0; iFF<dat::numsfdata; ++iFF)
	{
		std::ostringstream ostr;
		ostr << "ffacts.atom_" << iFF;
		std::string strAtom = ostr.str();

		std::string strA, strB;
		for(int i=0; i<5; ++i)
		{
			strA += tl::var_to_str(dat::sfdata[iFF].a[i], g_iPrec) + " ";
			strB += tl::var_to_str(dat::sfdata[iFF].b[i], g_iPrec) + " ";
		}

		prop.Add(strAtom + ".name", std::string(dat::sfdata[iFF].atomname));
		prop.Add(strAtom + ".a", strA);
		prop.Add(strAtom + ".b", strB);
		prop.Add(strAtom + ".c", tl::var_to_str(dat::sfdata[iFF].c, g_iPrec));
	}


	if(!prop.Save("res/data/ffacts.xml.gz"))
	{
		tl::log_err("Cannot write \"res/data/ffacts.xml.gz\".");
		return false;
	}
	return true;
}


// ============================================================================


bool gen_spacegroups_clp()
{
	tl::Prop<std::string> prop;
	prop.SetSeparator('.');

	const unsigned int iNumSGs = 230;
	prop.Add("sgroups.source", "Space group data extracted from Clipper (by K. Cowtan).");
	prop.Add("sgroups.source_url", "http://www.ysbl.york.ac.uk/~cowtan/clipper/");
	prop.Add("sgroups.num_groups", tl::var_to_str(iNumSGs));

	for(unsigned int iSG=1; iSG<=iNumSGs; ++iSG)
	{
		xtl::SpaceGroupClp sg(iSG);

		std::ostringstream ostr;
		ostr << "sgroups.group_" << (iSG-1);
		std::string strGroup = ostr.str();

		prop.Add(strGroup + ".number", tl::var_to_str(iSG));
		prop.Add(strGroup + ".name", sg.GetName());
		//prop.Add(strGroup + ".pointgroup", get_pointgroup(sg.GetName()));
		prop.Add(strGroup + ".lauegroup", sg.GetLaueGroup());
		//prop.Add(strGroup + ".crystalsys", sg.GetCrystalSystem());
		//prop.Add(strGroup + ".crystalsysname", sg.GetCrystalSystemName());


		std::vector<t_mat> vecTrafos, vecInv, vecPrim, vecCenter;
		sg.GetSymTrafos(vecTrafos);
		sg.GetInvertingSymTrafos(vecInv);
		sg.GetPrimitiveSymTrafos(vecPrim);
		sg.GetCenteringSymTrafos(vecCenter);


		prop.Add(strGroup + ".num_trafos", tl::var_to_str(vecTrafos.size()));
		unsigned int iTrafo = 0;
		for(const t_mat& matTrafo : vecTrafos)
		{
			bool bIsInv = xtl::is_mat_in_container<std::vector, t_mat>(vecInv, matTrafo);
			bool bIsPrim = xtl::is_mat_in_container<std::vector, t_mat>(vecPrim, matTrafo);
			bool bIsCenter = xtl::is_mat_in_container<std::vector, t_mat>(vecCenter, matTrafo);

			std::string strOpts = "; ";
			if(bIsPrim) strOpts += "p";
			if(bIsInv) strOpts += "i";
			if(bIsCenter) strOpts += "c";

			std::ostringstream ostrTrafo;
			ostrTrafo << strGroup << ".trafo_" << iTrafo;
			std::string strTrafo = ostrTrafo.str();

			prop.Add(strTrafo, tl::var_to_str(matTrafo, g_iPrec) + strOpts);

			++iTrafo;
		}
	}


	if(!prop.Save("res/data/sgroups.xml.gz"))
	{
		tl::log_err("Cannot write \"res/data/sgroups.xml.gz\".");
		return false;
	}

	return true;
}




// ============================================================================
