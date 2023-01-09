/**
 * S(q,w) module example for a pre-calculated grid (file format version 2)
 * @author Tobias Weber <tweber@ill.fr>
 * @date 06-jan-2020
 * @license GPLv2
 */

// Module: g++ -std=c++11 -shared -fPIC -o plugins/sqwgrid.so -I. -I/usr/include/QtCore -I/usr/include/x86_64-linux-gnu/qt5 -I/usr/include/x86_64-linux-gnu/qt5/QtCore examples/sqw_grid/grid_ver2.cpp tools/monteconvo/sqwbase.cpp tlibs/log/log.cpp tlibs/string/eval.cpp
// Test: g++ -std=c++11 -DDO_TEST -fPIC -o sqwgrid-tst -I. -I/usr/include/QtCore -I/usr/include/x86_64-linux-gnu/qt5 -I/usr/include/x86_64-linux-gnu/qt5/QtCore examples/sqw_grid/grid_ver2.cpp tools/monteconvo/sqwbase.cpp tlibs/string/eval.cpp tlibs/log/log.cpp -lboost_iostreams -lQt5Core

#include "grid_ver2.h"

#include "libs/version.h"
#include "tlibs/string/string.h"
#include "tlibs/math/math.h"
#include "tlibs/phys/neutrons.h"

#include <QFile>


using t_real = typename SqwMod::t_real;


// ----------------------------------------------------------------------------
// constructors

SqwMod::SqwMod()
{
	SqwBase::m_bOk = 0;
}


SqwMod::SqwMod(const std::string& strDatFile) : m_strDataFile(strDatFile)
{
	SqwBase::m_bOk = 0;
	tl::log_info("Loading grid version 2 data file: \"", strDatFile, "\".");


	QFile fileIdx(strDatFile.c_str());
	if(!fileIdx.exists())
	{
		tl::log_err("Grid data file \"", strDatFile, "\" does not exist.");
		return;
	}

	if(!fileIdx.open(QIODevice::ReadOnly))
	{
		tl::log_err("Grid data file \"", strDatFile, "\" cannot be opened.");
		return;
	}

	const void *pMemIdx = fileIdx.map(0, sizeof(m_indexBlockOffset) + 9*sizeof(t_real));
	if(!pMemIdx)
	{
		tl::log_err("Grid data file \"", strDatFile, "\" cannot be mapped.");
		return;
	}

	m_indexBlockOffset = *((std::size_t*)pMemIdx);

	m_hmin = ((t_real*)pMemIdx)[1];
	m_hmax = ((t_real*)pMemIdx)[2];
	m_hstep = ((t_real*)pMemIdx)[3];
	m_hsize = std::size_t(std::round((m_hmax - m_hmin) / m_hstep));

	m_kmin = ((t_real*)pMemIdx)[4];
	m_kmax = ((t_real*)pMemIdx)[5];
	m_kstep = ((t_real*)pMemIdx)[6];
	m_ksize = std::size_t(std::round((m_kmax - m_kmin) / m_kstep));

	m_lmin = ((t_real*)pMemIdx)[7];
	m_lmax = ((t_real*)pMemIdx)[8];
	m_lstep = ((t_real*)pMemIdx)[9];
	m_lsize = std::size_t(std::round((m_lmax - m_lmin) / m_lstep));

	std::size_t numEntries =
		std::size_t(((m_hmax-m_hmin) / m_hstep)) *
		std::size_t(((m_kmax-m_kmin) / m_kstep)) *
		std::size_t(((m_lmax-m_lmin) / m_lstep));

	fileIdx.unmap((unsigned char*)pMemIdx);
	fileIdx.close();

	tl::log_info("Data block dimensions: h=", m_hmin, "..", m_hmax, " (delta=", m_hstep, "), ",
		"k=", m_kmin, "..", m_kmax, " (delta=", m_kstep, "), ",
		"l=", m_lmin, "..", m_lmax, " (delta=", m_lstep, ").");
	tl::log_info("Number of entries: ", numEntries, ".");
	//tl::log_info("Offset of index block: ", m_indexBlockOffset, ".");

	SqwBase::m_bOk = 1;
}


SqwMod::~SqwMod()
{
}


// ----------------------------------------------------------------------------
// dispersion, spectral weight and structure factor

std::tuple<std::vector<t_real>, std::vector<t_real>>
	SqwMod::disp(t_real dh, t_real dk, t_real dl) const
{
	/**
	 * calculate file index based on coordinates
	 */
	auto hkl_to_idx = [this](t_real h, t_real k, t_real l) -> std::size_t
	{
		// clamp values to boundaries
		if(h < m_hmin) h = m_hmin;
		if(k < m_kmin) k = m_kmin;
		if(l < m_lmin) l = m_lmin;
		if(h >= m_hmax) h = m_hmax - m_hstep;
		if(k >= m_kmax) k = m_kmax - m_kstep;
		if(l >= m_lmax) l = m_lmax - m_lstep;

		// position indices
		std::size_t iH = std::size_t(std::round(((h - m_hmin) / m_hstep)));
		std::size_t iK = std::size_t(std::round(((k - m_kmin) / m_kstep)));
		std::size_t iL = std::size_t(std::round(((l - m_lmin) / m_lstep)));

		// clamp again
		if(iH >= m_hsize) iH = m_hsize-1;
		if(iK >= m_ksize) iK = m_ksize-1;
		if(iL >= m_lsize) iL = m_lsize-1;

		return iH*m_ksize*m_lsize + iK*m_lsize + iL;
	};



	QFile fileDat(m_strDataFile.c_str());
	if(!fileDat.exists())
	{
		tl::log_err("Grid data file \"", m_strDataFile, "\" does not exist.");
		return std::make_tuple(std::vector<t_real>(), std::vector<t_real>());
	}

	if(!fileDat.open(QIODevice::ReadOnly))
	{
		tl::log_err("Data file \"", m_strDataFile, "\" cannot be opened.");
		return std::make_tuple(std::vector<t_real>(), std::vector<t_real>());
	}


	// ------------------------------------------------------------------------
	// the index block the offsets into the data block
	std::size_t idx_file_offs = hkl_to_idx(dh, dk, dl);

	const void *pMemIdx = fileDat.map(
		m_indexBlockOffset + idx_file_offs*sizeof(std::size_t), sizeof(std::size_t));
	if(!pMemIdx)
	{
		tl::log_err("Grid data file \"", m_strDataFile, "\" cannot be mapped.");
		return std::make_tuple(std::vector<t_real>(), std::vector<t_real>());
	}

	std::size_t dat_file_offs = *((std::size_t*)pMemIdx);
	fileDat.unmap((unsigned char*)pMemIdx);
	// ------------------------------------------------------------------------


	// ------------------------------------------------------------------------
	// the data block holds the energies and spectral weights of the dispersion branches

	const void *pMemDat = fileDat.map(dat_file_offs, sizeof(std::size_t));
	if(!pMemDat)
	{
		tl::log_err("Grid data file \"", m_strDataFile, "\" cannot be mapped (1).");
		return std::make_tuple(std::vector<t_real>(), std::vector<t_real>());
	}

	// number of dispersion branches and weights
	unsigned int iNumBranches = *((unsigned int*)pMemDat);
	fileDat.unmap((unsigned char*)pMemDat);


	// map actual (E, w) data
	const t_real *pBranches = (t_real*)fileDat.map(dat_file_offs+sizeof(iNumBranches),
		iNumBranches*sizeof(t_real)*2);
	if(!pMemDat)
	{
		tl::log_err("Grid data file \"", m_strDataFile, "\" cannot be mapped (2).");
		return std::make_tuple(std::vector<t_real>(), std::vector<t_real>());
	}


	std::vector<t_real> vecE, vecw;
	for(unsigned int iBranch=0; iBranch<iNumBranches; ++iBranch)
	{
		if(!tl::float_equal(pBranches[iBranch*2 + 1], t_real(0)))
		{
			vecE.push_back(pBranches[iBranch*2 + 0]);	// energy
			vecw.push_back(pBranches[iBranch*2 + 1]);	// weight
		}
	}

	fileDat.unmap((unsigned char*)pBranches);
	fileDat.close();
	// ------------------------------------------------------------------------


	return std::make_tuple(vecE, vecw);
}


/**
 * S(q,E)
 */
t_real SqwMod::operator()(t_real dh, t_real dk, t_real dl, t_real dE) const
{
	std::vector<t_real> vecE, vecW;
	std::tie(vecE, vecW) = disp(dh, dk, dl);

	t_real dInc=0, dS_p=0, dS_m=0;
	if(!tl::float_equal(m_dIncAmp, t_real(0)))
		dInc = tl::gauss_model(dE, t_real(0), m_dIncSigma, m_dIncAmp, t_real(0));

	t_real dS = 0;
	for(std::size_t iE=0; iE<vecE.size(); ++iE)
		dS += tl::gauss_model(dE, vecE[iE], m_dSigma, vecW[iE], t_real(0));

	return m_dS0*dS * tl::bose_cutoff(dE, m_dT, m_dcut) + dInc;
}



// ----------------------------------------------------------------------------
// get & set variables

std::vector<SqwMod::t_var> SqwMod::GetVars() const
{
	std::vector<t_var> vecVars;

	vecVars.push_back(SqwBase::t_var{"T", "real", tl::var_to_str(m_dT)});
	vecVars.push_back(SqwBase::t_var{"bose_cutoff", "real", tl::var_to_str(m_dcut)});
	vecVars.push_back(SqwBase::t_var{"sigma", "real", tl::var_to_str(m_dSigma)});
	vecVars.push_back(SqwBase::t_var{"inc_amp", "real", tl::var_to_str(m_dIncAmp)});
	vecVars.push_back(SqwBase::t_var{"inc_sigma", "real", tl::var_to_str(m_dIncSigma)});
	vecVars.push_back(SqwBase::t_var{"S0", "real", tl::var_to_str(m_dS0)});

	return vecVars;
}


void SqwMod::SetVars(const std::vector<SqwMod::t_var>& vecVars)
{
	if(!vecVars.size()) return;

	for(const SqwBase::t_var& var : vecVars)
	{
		const std::string& strVar = std::get<0>(var);
		const std::string& strVal = std::get<2>(var);

		if(strVar == "T") m_dT = tl::str_to_var<decltype(m_dT)>(strVal);
		else if(strVar == "bose_cutoff") m_dcut = tl::str_to_var<decltype(m_dcut)>(strVal);
		else if(strVar == "sigma") m_dSigma = tl::str_to_var<t_real>(strVal);
		else if(strVar == "inc_amp") m_dIncAmp = tl::str_to_var<decltype(m_dIncAmp)>(strVal);
		else if(strVar == "inc_sigma") m_dIncSigma = tl::str_to_var<decltype(m_dIncSigma)>(strVal);
		else if(strVar == "S0") m_dS0 = tl::str_to_var<decltype(m_dS0)>(strVal);
	}
}


bool SqwMod::SetVarIfAvail(const std::string& strKey, const std::string& strNewVal)
{
	return SqwBase::SetVarIfAvail(strKey, strNewVal);
}



// ----------------------------------------------------------------------------
// copy

SqwBase* SqwMod::shallow_copy() const
{
	SqwMod *pMod = new SqwMod();

	pMod->m_dT = this->m_dT;
	pMod->m_dcut = this->m_dcut;
	pMod->m_dSigma = this->m_dSigma;
	pMod->m_dIncAmp = this->m_dIncAmp;
	pMod->m_dIncSigma = this->m_dIncSigma;
	pMod->m_dS0 = this->m_dS0;

	pMod->m_strDataFile = this->m_strDataFile;
	pMod->m_indexBlockOffset = this->m_indexBlockOffset;

	pMod->m_hmin = this->m_hmin;
	pMod->m_hmax = this->m_hmax;
	pMod->m_hstep = this->m_hstep;
	pMod->m_hsize = this->m_hsize;

	pMod->m_kmin = this->m_kmin;
	pMod->m_kmax = this->m_kmax;
	pMod->m_kstep = this->m_kstep;
	pMod->m_ksize = this->m_ksize;

	pMod->m_lmin = this->m_lmin;
	pMod->m_lmax = this->m_lmax;
	pMod->m_lstep = this->m_lstep;
	pMod->m_lsize = this->m_lsize;

	return pMod;
}



#ifdef DO_TEST
// ----------------------------------------------------------------------------
// test querying of data

int main(int argc, char **argv)
{
	if(argc <= 1)
	{
		std::cerr << "Please specify a grid data file." << std::endl;
		return -1;
	}


	SqwMod mod(argv[1]);
	//mod.SetVarIfAvail("sigma", "0.1");

	while(1)
	{
		t_real h, k, l;
		std::cout << "Enter hkl: ";
		std::cin >> h >> k >> l;

		std::vector<t_real> vecE, vecW;
		std::tie(vecE, vecW) = mod.disp(h, k, l);

		for(std::size_t i=0; i<vecE.size(); ++i)
			std::cout << "(" << i+1 << ") " << "E = " << vecE[i] << ", weight = " << vecW[i] << "\n";
		std::cout << std::endl;
	}
}
// ----------------------------------------------------------------------------

#else

// ----------------------------------------------------------------------------
// SO interface
#include <boost/dll/alias.hpp>


static const char* pcModIdent = "grid2mod";
static const char* pcModName = "Grid Version 2";


std::tuple<std::string, std::string, std::string> sqw_info()
{
	//tl::log_info("In ", __func__, ".");
	return std::make_tuple(TAKIN_VER, pcModIdent, pcModName);
}


std::shared_ptr<SqwBase> sqw_construct(const std::string& strCfgFile)
{
	//tl::log_info("In ", __func__, ".");
	return std::make_shared<SqwMod>(strCfgFile);
}


// exports from so file
BOOST_DLL_ALIAS(sqw_info, takin_sqw_info);
BOOST_DLL_ALIAS(sqw_construct, takin_sqw);
// ----------------------------------------------------------------------------

#endif
