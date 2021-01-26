// gcc -O2 -march=native -DNDEBUG -o export export.cpp -std=c++11 -lstdc++ -lm -I../.. ../../libs/spacegroups/spacegroup.cpp ../../libs/spacegroups/crystalsys.cpp ../../libs/formfactors/formfact.cpp ../../tlibs/log/log.cpp ../../libs/globals.cpp -DNO_QT -lboost_system -lboost_filesystem -lboost_iostreams
/**
 * exports tables
 * @author Tobias Weber <tobias.weber@tum.de>
 * @date nov-2017
 * @license GPLv2
 */

#include "libs/spacegroups/spacegroup.h"
#include "libs/formfactors/formfact.h"

namespace ublas = tl::ublas;
using t_real = double;
std::size_t g_iprec = std::numeric_limits<t_real>::max_digits10 - 2;


static void export_jl(std::shared_ptr<const xtl::SpaceGroups<t_real>> sgs,
	std::shared_ptr<const xtl::ScatlenList<t_real>> lst,
	std::shared_ptr<const xtl::FormfactList<t_real>> lstff,
	std::shared_ptr<const xtl::MagFormfactList<t_real>> lstmff)
{
	// Names of identifiers
	const char* pcNum = "num";
	const char* pcName = "name";
	const char* pcTrafo = "trafos";
	const char* pcCoh = "coh";
	const char* pcInc = "inc";

	// --------------------------------------------------------------------------------------------
	// Space groups
	// --------------------------------------------------------------------------------------------
	tl::log_info("Writing space groups to jl ...");
	std::ofstream ofstrSGs("sgs.jl");
	ofstrSGs.precision(g_iprec);

	ofstrSGs << "sgs = Dict(\n";

	const auto* pSGs = sgs->get_space_groups_vec();
	for(const auto* pSg : *pSGs)
	{
		ofstrSGs << "\"" << pSg->GetName() << "\" => Dict(";
		ofstrSGs << "\"" << pcNum << "\" => " << pSg->GetNr();
		ofstrSGs << ", \"" << pcName << "\" => \"" << pSg->GetName() << "\"";
		ofstrSGs << ", \"" << pcTrafo << "\" => \n";
		ofstrSGs << "[\n";

		const auto& vecTrafos = pSg->GetTrafos();
		for(const auto& mat : vecTrafos)
		{
			ofstrSGs << "\t[ "
				<< mat(0,0) << " " << mat(0,1) << " " << mat(0,2) << " " << mat(0,3) << "; "
				<< mat(1,0) << " " << mat(1,1) << " " << mat(1,2) << " " << mat(1,3) << "; "
				<< mat(2,0) << " " << mat(2,1) << " " << mat(2,2) << " " << mat(2,3) << "; "
				<< mat(3,0) << " " << mat(3,1) << " " << mat(3,2) << " " << mat(3,3)
				<< " ],\n";
		}

		ofstrSGs << "]";
		ofstrSGs << "),\n";
	}

	ofstrSGs << ");\n";

	// Test
	//ofstrSGs << "\n\nprintln(sgs[\"P2_13\"]);\n";
	// --------------------------------------------------------------------------------------------


	// --------------------------------------------------------------------------------------------
	// Scattering lengths
	// --------------------------------------------------------------------------------------------
	tl::log_info("Writing scattering lengths to jl ...");
	std::ofstream ofstrSlens("slens.jl");
	ofstrSlens.precision(g_iprec);

	ofstrSlens << "slens = Dict(\n";

	using sl_type = typename xtl::ScatlenList<t_real>::elem_type;
	auto write_elem = [pcName, pcCoh, pcInc]
	(std::ostream& ostr, const sl_type& elem)
	{
		ostr << "\"" << elem.GetAtomIdent() << "\" => Dict( ";
		ostr << "\"" << pcName << "\" => \"" << elem.GetAtomIdent() << "\"";
		ostr << ", \"" << pcCoh << "\" => " << elem.GetCoherent().real() <<
		" + " << elem.GetCoherent().imag() << "im";
		ostr << ", \"" << pcInc << "\" => " << elem.GetIncoherent().real() <<
		" + " << elem.GetIncoherent().imag() << "im";
		ostr << " ),\n";
	};

	for(std::size_t iElem=0; iElem<lst->GetNumElems(); ++iElem)
		write_elem(ofstrSlens, lst->GetElem(iElem));
	for(std::size_t iElem=0; iElem<lst->GetNumIsotopes(); ++iElem)
		write_elem(ofstrSlens, lst->GetIsotope(iElem));

	ofstrSlens << ");\n";

	// Test
	//ofstrSlens << "\n\nprintln(slens[\"Mn\"][\"coh\"]);\n";
	// --------------------------------------------------------------------------------------------
}


static void export_js(std::shared_ptr<const xtl::SpaceGroups<t_real>> sgs,
	std::shared_ptr<const xtl::ScatlenList<t_real>> lst,
	std::shared_ptr<const xtl::FormfactList<t_real>> lstff,
	std::shared_ptr<const xtl::MagFormfactList<t_real>> lstmff)
{
	// Names of identifiers
	const char* pcNum = "num";
	const char* pcName = "name";
	const char* pcTrafo = "trafos";
	const char* pcCohRe = "coh_r";
	const char* pcCohIm = "coh_i";
	const char* pcIncRe = "inc_r";
	const char* pcIncIm = "inc_i";

	// --------------------------------------------------------------------------------------------
	// Space groups
	// --------------------------------------------------------------------------------------------
	tl::log_info("Writing space groups to js ...");
	std::ofstream ofstrSGs("sgs.js");
	ofstrSGs.precision(g_iprec);

	ofstrSGs << "sgs = {\n";

	const auto* pSGs = sgs->get_space_groups_vec();
	for(const auto* pSg : *pSGs)
	{
		ofstrSGs << "\"" << pSg->GetName() << "\" : { ";
		ofstrSGs << pcNum << ": " << pSg->GetNr();
		ofstrSGs << ", " << pcName << ": \"" << pSg->GetName() << "\"";
		ofstrSGs << ", " << pcTrafo << ": \n";
		ofstrSGs << "[\n";

		const auto& vecTrafos = pSg->GetTrafos();
		for(const auto& mat : vecTrafos)
		{
			ofstrSGs << "\t[ "
				<< mat(0,0) << "," << mat(0,1) << "," << mat(0,2) << "," << mat(0,3) << ", "
				<< mat(1,0) << "," << mat(1,1) << "," << mat(1,2) << "," << mat(1,3) << ", "
				<< mat(2,0) << "," << mat(2,1) << "," << mat(2,2) << "," << mat(2,3) << ", "
				<< mat(3,0) << "," << mat(3,1) << "," << mat(3,2) << "," << mat(3,3)
				<< " ],\n";
		}

		ofstrSGs << "]";
		ofstrSGs << "},\n";
	}

	ofstrSGs << "};\n";

	// Test
	//ofstrSGs << "\n\nconsole.log(sgs[\"P2_13\"]);\n";
	// --------------------------------------------------------------------------------------------


	// --------------------------------------------------------------------------------------------
	// Scattering lengths
	// --------------------------------------------------------------------------------------------
	tl::log_info("Writing scattering lengths to js ...");
	std::ofstream ofstrSlens("slens.js");
	ofstrSlens.precision(g_iprec);

	ofstrSlens << "slens = {\n";

	using sl_type = typename xtl::ScatlenList<t_real>::elem_type;
	auto write_elem = [pcName, pcCohRe, pcCohIm, pcIncRe, pcIncIm]
	(std::ostream& ostr, const sl_type& elem)
	{
		ostr << "\"" << elem.GetAtomIdent() << "\" : { ";
		ostr << pcName << ": \"" << elem.GetAtomIdent() << "\"";
		ostr << ", " << pcCohRe << ": " << elem.GetCoherent().real();
		ostr << ", " << pcCohIm << ": " << elem.GetCoherent().imag();
		ostr << ", " << pcIncRe << ": " << elem.GetIncoherent().real();
		ostr << ", " << pcIncIm << ": " << elem.GetIncoherent().imag();
		ostr << " },\n";
	};

	for(std::size_t iElem=0; iElem<lst->GetNumElems(); ++iElem)
		write_elem(ofstrSlens, lst->GetElem(iElem));
	for(std::size_t iElem=0; iElem<lst->GetNumIsotopes(); ++iElem)
		write_elem(ofstrSlens, lst->GetIsotope(iElem));

	ofstrSlens << "};\n";

	// Test
	//ofstrSlens << "\n\nconsole.log(slens[\"Mn\"].coh_r);\n";
	// --------------------------------------------------------------------------------------------
}


int main(int argc, char** argv)
{
	// --------------------------------------------------------------------------------------------
	// Load tables
	// --------------------------------------------------------------------------------------------
	std::shared_ptr<const xtl::SpaceGroups<t_real>> sgs = xtl::SpaceGroups<t_real>::GetInstance();
	std::shared_ptr<const xtl::ScatlenList<t_real>> lst = xtl::ScatlenList<t_real>::GetInstance();
	std::shared_ptr<const xtl::FormfactList<t_real>> lstff = xtl::FormfactList<t_real>::GetInstance();
	std::shared_ptr<const xtl::MagFormfactList<t_real>> lstmff = xtl::MagFormfactList<t_real>::GetInstance();
	// --------------------------------------------------------------------------------------------

	export_jl(sgs, lst, lstff, lstmff);
	export_js(sgs, lst, lstff, lstmff);
}
