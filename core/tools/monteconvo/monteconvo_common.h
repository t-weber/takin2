/**
 * monte carlo convolution tool -- common functions
 * @author Tobias Weber <tweber@ill.fr>
 * @date sep-2020
 * @license GPLv2
 */

#ifndef __MONTECONVO_COMMON_H__
#define __MONTECONVO_COMMON_H__

#include "tlibs/math/math.h"
#include "tlibs/string/string.h"

#include "tools/convofit/scan.h"
#include "TASReso.h"


#define EPS_RLU 1e-3



/**
 * determines the x axis of the scan
 */
template<class t_real = double>
std::tuple<bool, int, std::string, std::vector<std::vector<t_real>>>
get_scan_axis(bool bIncludeE, int iScanAxis, std::size_t num_steps, t_real eps,
    t_real startH, t_real endH, t_real startK, t_real endK, t_real startL, t_real endL,
    t_real startE, t_real endE)
{
	std::vector<t_real> vecH = tl::linspace<t_real,t_real>(startH, endH, num_steps);
	std::vector<t_real> vecK = tl::linspace<t_real,t_real>(startK, endK, num_steps);
	std::vector<t_real> vecL = tl::linspace<t_real,t_real>(startL, endL, num_steps);
	std::vector<t_real> vecE = tl::linspace<t_real,t_real>(startE, endE, num_steps);

	std::vector<std::vector<t_real>> vecScanAxes{{ std::move(vecH), std::move(vecK), std::move(vecL) }};
	if(bIncludeE)
		vecScanAxes.emplace_back(std::move(vecE));

	int iScanAxisIdx = 0;

	std::string strScanVar = "";
	// either scan axis is directly selected OR (automatic is set AND the start/stop values are different)
	if(iScanAxis==1 || (iScanAxis==0 && !tl::float_equal(startH, endH, eps)))
	{
		strScanVar = "h (rlu)";
		iScanAxisIdx = 0;
	}
	else if(iScanAxis==2 || (iScanAxis==0 && !tl::float_equal(startK, endK, eps)))
	{
		strScanVar = "k (rlu)";
		iScanAxisIdx = 1;
	}
	else if(iScanAxis==3 || (iScanAxis==0 && !tl::float_equal(startL, endL, eps)))
	{
		strScanVar = "l (rlu)";
		iScanAxisIdx = 2;
	}
	else if(bIncludeE && (iScanAxis==4 || (iScanAxis==0 && !tl::float_equal(startE, endE, eps))))
	{
		strScanVar = "E (meV)";
		iScanAxisIdx = 3;
	}
	else
	{
		return std::make_tuple(false, iScanAxisIdx, strScanVar, vecScanAxes);
	}

	return std::make_tuple(true, iScanAxisIdx, strScanVar, vecScanAxes);
}



extern ResoFocus get_reso_focus(int iFocMono, int iFocAna);


extern bool load_scan_file(const std::string& _strFile, Scan& scan,
	bool bFlipAxis=false, const Filter& filter=Filter());


#endif
