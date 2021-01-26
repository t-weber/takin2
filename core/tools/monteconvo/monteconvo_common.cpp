/**
 * monte carlo convolution tool -- common functions
 * @author Tobias Weber <tweber@ill.fr>
 * @date sep-2020
 * @license GPLv2
 */

#include "monteconvo_common.h"


#define EPS_RLU 1e-3


ResoFocus get_reso_focus(int iFocMono, int iFocAna)
{
	unsigned ifocMode = unsigned(ResoFocus::FOC_UNCHANGED);
	if(iFocMono == 1) ifocMode |= unsigned(ResoFocus::FOC_MONO_FLAT);	// flat
	else if(iFocMono == 2) ifocMode |= unsigned(ResoFocus::FOC_MONO_H);	// horizontal
	else if(iFocMono == 3) ifocMode |= unsigned(ResoFocus::FOC_MONO_V);	// vertical
	else if(iFocMono == 4) ifocMode |= unsigned(ResoFocus::FOC_MONO_V)|unsigned(ResoFocus::FOC_MONO_H);		// both
	if(iFocAna == 1) ifocMode |= unsigned(ResoFocus::FOC_ANA_FLAT);		// flat
	else if(iFocAna == 2) ifocMode |= unsigned(ResoFocus::FOC_ANA_H);	// horizontal
	else if(iFocAna == 3) ifocMode |= unsigned(ResoFocus::FOC_ANA_V);	// vertical
	else if(iFocAna == 4) ifocMode |= unsigned(ResoFocus::FOC_ANA_V)|unsigned(ResoFocus::FOC_ANA_H);		// both

	return ResoFocus(ifocMode);
}



bool load_scan_file(const std::string& _strFile, Scan& scan,
	bool bFlipAxis, const Filter& filter)
{
	std::string strFile = _strFile;
	tl::trim(strFile);
	if(strFile == "")
		return false;

	std::vector<std::string> vecFiles;
	tl::get_tokens<std::string, std::string>(strFile, ";", vecFiles);
	std::for_each(vecFiles.begin(), vecFiles.end(), [](std::string& str){ tl::trim(str); });

	scan = Scan();
	bool bLoaded = ::load_file(vecFiles, scan, 1, filter, bFlipAxis);

	// if file was not found, alternatively look in global paths
	if(!bLoaded)
	{
		const std::vector<std::string>& vecGlobPaths = get_global_paths();
		for(const std::string& strGlobPath : vecGlobPaths)
		{
			std::vector<std::string> _vecFiles;
			for(const std::string& _strFile : vecFiles)
				_vecFiles.push_back(strGlobPath + "/" + _strFile);

			if((bLoaded = ::load_file(_vecFiles, scan, 1, filter, bFlipAxis)))
				break;
		}
	}

    return bLoaded;
}
