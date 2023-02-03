/**
 * invocation of gnuplot
 * @author Tobias Weber <tobias.weber@tum.de>
 * @date 24-dec-2013
 * @license GPLv2 or GPLv3
 *
 * ----------------------------------------------------------------------------
 * tlibs -- a physical-mathematical C++ template library
 * Copyright (C) 2017-2021  Tobias WEBER (Institut Laue-Langevin (ILL),
 *                          Grenoble, France).
 * Copyright (C) 2015-2017  Tobias WEBER (Technische Universitaet Muenchen
 *                          (TUM), Garching, Germany).
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 2 of the License, or
 * (at your option) version 3.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 * ----------------------------------------------------------------------------
 */

#ifndef __GPL_PLOTTER_IMPL_H__
#define __GPL_PLOTTER_IMPL_H__

#include "../helper/misc.h"
#include "../string/string.h"
#include "../log/log.h"

#include <cstdio>
#include <sstream>
#include <fstream>
#include <iomanip>

#include "gnuplot.h"


namespace tl {

template<class t_real>
void GnuPlot<t_real>::DeInit()
{
	m_iStartCounter = 0;
	m_pProc.reset(nullptr);
}


template<class t_real>
bool GnuPlot<t_real>::Init()
{
	if(IsReady()) return true;
	DeInit();
	const std::string strGplTool = "gnuplot";


	// get gnuplot version
	{
		PipeProc<char> proc((strGplTool + " 2>/dev/null --version").c_str(), 0);
		if(!proc.IsReady())
			return "";

		std::getline(proc.GetIstr(), m_strVersion);
		trim(m_strVersion);

		if(m_strVersion == "")
			return false;
	}


	// open command pipe to gnuplot
	{
		m_pProc.reset(new PipeProc<char>((strGplTool + " -p 2>/dev/null 1>/dev/null").c_str(), 1));
		if(!m_pProc || !m_pProc->IsReady())
			return false;

		m_pProc->GetOstr().precision(m_iPrec);
		m_pProc->GetOstr() << "set grid\n";
		m_pProc->GetOstr() << "set nokey\n";
		m_pProc->GetOstr() << "set size 1,1\n";
		m_pProc->GetOstr() << "set palette rgbformulae 33,13,10\n";
	}

	return true;
}


template<class t_real>
void GnuPlot<t_real>::SetTerminal(int iWnd, const char* pcBackend, t_real dW, t_real dH)
{
	if(!IsReady()) return;
	if(m_bTermLocked) return;

	m_pProc->GetOstr() << "set output\n";
	m_pProc->GetOstr() << "set obj 1 rectangle behind fillcolor rgb \"#ffffff\" from screen 0,0 to screen 1,1\n";

	m_pProc->GetOstr() << "set term " << pcBackend <<  " " << iWnd << " "
		<< "enhanced font \"NimbusSanL-Regu,12\" persist dashed ";

	if(dW>=t_real(0) && dH>=t_real(0))
		m_pProc->GetOstr() << "size " << dW << "," << dH << " ";
	else
		m_pProc->GetOstr() << "size 640,480 ";

	m_pProc->GetOstr() << "\n";
}


template<class t_real>
void GnuPlot<t_real>::SetFileTerminal(const char* pcFile, t_real dW, t_real dH)
{
	if(!IsReady()) return;
	if(m_bTermLocked) return;

	std::string strFile = pcFile;
	std::string strExt = get_fileext(strFile);

	if(str_is_equal(strExt, std::string("pdf"), 0))
	{
		m_pProc->GetOstr() << "set term pdf "
			<< "enhanced color font \"NimbusSanL-Regu,16\" ";

		if(dW>=t_real(0) && dH>=t_real(0))
			m_pProc->GetOstr() << "size " << dW << "," << dH << " ";
		m_pProc->GetOstr() << "\n";
	}
	else if(str_is_equal(strExt, std::string("ps"), 0))
	{
		m_pProc->GetOstr() << "set term postscript eps "
			<< "enhanced color font \"NimbusSanL-Regu,16\" ";

		if(dW>=t_real(0) && dH>=t_real(0))
			m_pProc->GetOstr() << "size " << dW << "," << dH << " ";
		m_pProc->GetOstr() << "\n";
	}
	else
	{
		log_err("Unknown file extension \"", strExt, "\" for output terminal.");
		return;
	}

	m_pProc->GetOstr() << "set output \"" << strFile << "\"\n";
}


template<class t_real>
void GnuPlot<t_real>::SetPrec(unsigned int iPrec)
{
	if(!IsReady()) return;

	m_iPrec = iPrec;
	if(m_pProc)
		m_pProc->GetOstr().precision(m_iPrec);
}


template<class t_real>
void GnuPlot<t_real>::RefreshVars()
{
	if(!IsReady()) return;

	if(m_bHasLegend)
		m_pProc->GetOstr() << "set key on " << m_strLegendPlacement
			<< " box 1 " << m_strLegendOpts << "\n";
	else
		m_pProc->GetOstr() << "set nokey\n";
}


template<class t_real>
void GnuPlot<t_real>::SimplePlot(const std::vector<t_real>& vecX, const std::vector<t_real>& vecY,
	const std::vector<t_real>& vecYErr, const std::vector<t_real>& vecXErr,
	LineStyle style)
{
	if(!IsReady()) return;
	m_pProc->GetOstr() << "plot \"-\" ";

	switch(style)
	{
		case STYLE_LINES_SOLID:
			m_pProc->GetOstr() << "using ($1):($2) with lines linetype 1 linewidth 1 ";
			break;
		case STYLE_LINES_DASHED:
			m_pProc->GetOstr() << "using ($1):($2) with lines linetype 2 linewidth 1 ";
			break;
		default:
		case STYLE_POINTS:
			break;
	}
	m_pProc->GetOstr() << "\n";

	std::string strTable = BuildTable(vecX, vecY, vecYErr, vecXErr);
	m_pProc->GetOstr() << strTable;
	m_pProc->GetOstr().flush();
}


template<class t_real>
void GnuPlot<t_real>::SimplePlot2d(const std::vector<std::vector<t_real> >& vec,
	t_real dMinX, t_real dMaxX, t_real dMinY, t_real dMaxY)
{
	if(!IsReady()) return;

	std::vector<std::size_t> vecSizes;
	vecSizes.reserve(vec.size());

	for(std::size_t iY=0; iY<vec.size(); ++iY)
		vecSizes.push_back(vec[iY].size());

	std::vector<std::size_t>::iterator iterMin =
		std::min_element(vecSizes.begin(), vecSizes.end());
	std::size_t iXCntMin = *iterMin;


	std::size_t iYDim = vec.size();
	std::size_t iXDim = iXCntMin;

	// invalid values select image dimensions
	if(dMinX > dMaxX)
	{
		dMinX = 0.;
		dMaxX = iXDim-1;
	}
	if(dMinY > dMaxY)
	{
		dMinY = 0.;
		dMaxY = iYDim-1;
	}

	// ----------------------------------------
	// ranges
	m_pProc->GetOstr() << "set tics out scale 0.75\n";

	t_real dRangeMinX = tic_trafo<t_real>(iXDim, dMinX, dMaxX, 0, -0.5);
	t_real dRangeMaxX = tic_trafo<t_real>(iXDim, dMinX, dMaxX, 0, t_real(iXDim)-0.5);
	t_real dRangeMinY = tic_trafo<t_real>(iYDim, dMinY, dMaxY, 0, -0.5);
	t_real dRangeMaxY = tic_trafo<t_real>(iYDim, dMinY, dMaxY, 0, t_real(iYDim)-0.5);

	m_pProc->GetOstr() << "set xrange [" << dRangeMinX << ":" << dRangeMaxX << "]\n";
	m_pProc->GetOstr() << "set yrange [" << dRangeMinY << ":" << dRangeMaxY << "]\n";
	// ----------------------------------------

	// ----------------------------------------
	// tics
	std::ostringstream ostrTics;
	ostrTics.precision(m_iPrec);

	ostrTics << "using (" << dMinX << " + " << "($1)/" << iXDim
		<< " * (" << dMaxX << "-" << dMinX << "))" << " : "
		<< "(" << dMinY << " + " << "($2)/" << iYDim
		<< " * (" << dMaxY << "-" << dMinY << "))" << " : ($3)";

	std::string strTics = ostrTics.str();
	// ----------------------------------------

	m_pProc->GetOstr() << "plot \"-\" " << strTics << " matrix with image\n";


	for(std::size_t iY=0; iY<vec.size(); ++iY)
	{
		for(std::size_t iX=0; iX<iXCntMin; ++iX)
			m_pProc->GetOstr() << vec[iY][iX] << " ";
		m_pProc->GetOstr() << "\n";
	}

	m_pProc->GetOstr() << "end\nend\n";


	m_pProc->GetOstr().flush();
}


template<class t_real>
void GnuPlot<t_real>::AddLine(const PlotObj<t_real>& obj)
{
	m_vecObjs.push_back(obj);
}


template<class t_real>
void GnuPlot<t_real>::StartPlot()
{
	if(!IsReady()) return;

	if(m_iStartCounter == 0)
		m_vecObjs.clear();

	++m_iStartCounter;
}


template<class t_real>
void GnuPlot<t_real>::SetCmdFileOutput(const char* pcFile)
{
	m_strCmdFileOutput = pcFile;
}


template<class t_real>
void GnuPlot<t_real>::FinishPlot()
{
	if(!IsReady()) return;

	if(--m_iStartCounter == 0)
	{
		std::string strCmd = BuildCmd();
		RefreshVars();

		if(m_strCmdFileOutput != "")
		{
			std::ofstream ofCmd(m_strCmdFileOutput);
			if(!!ofCmd)
			{
				ofCmd << "#!/usr/bin/gnuplot -p\n\n";
				ofCmd << "#set term pdf enhanced color font \"NimbusSanL-Regu,16\"\n";
				ofCmd << "#set output \"" << m_strCmdFileOutput << ".pdf\"\n\n";

				ofCmd << strCmd << "\n";
				ofCmd.close();
				m_strCmdFileOutput = "";
			}
		}

		m_pProc->GetOstr() << strCmd;
		m_pProc->GetOstr().flush();
		m_vecObjs.clear();
	}
}


template<class t_real>
std::string GnuPlot<t_real>::BuildCmd()
{
	m_bHasLegend = 0;

	std::ostringstream ostr;
	ostr.precision(m_iPrec);

	ostr << "plot ";

	for(const PlotObj<t_real>& obj : m_vecObjs)
	{
		const bool bConnectLines = (obj.linestyle != STYLE_POINTS);
		const bool bHasXErr = (obj.vecErrX.size() != 0);
		const bool bHasYErr = (obj.vecErrY.size() != 0);
		t_real dSize = boost::get_optional_value_or(obj.odSize, bConnectLines ? 1. : 1.25);

		std::ostringstream ostrTmp;
		std::string strPointStyle;

		if(bHasXErr && bHasYErr)
			ostrTmp << "using ($1):($2):($3):($4) with xyerrorbars";
		else if(bHasXErr && !bHasYErr)
			ostrTmp << "using ($1):($2):($3) with xerrorbars";
		else if(!bHasXErr && bHasYErr)
			ostrTmp << "using ($1):($2):($3) with yerrorbars";
		else if(!bHasXErr && !bHasYErr)
			ostrTmp << "using ($1):($2) with points";

		ostrTmp << " pointtype 7 pointsize " << dSize;
		strPointStyle = ostrTmp.str();


		ostr << "\"-\" ";
		switch(obj.linestyle)
		{
			case STYLE_LINES_SOLID:
				ostr << "using ($1):($2) with lines linetype 1 linewidth " << dSize << " ";
				break;
			case STYLE_LINES_DASHED:
				ostr << "using ($1):($2) with lines linetype 2 linewidth " << dSize << " ";
				break;
			default:
				log_warn("Unknown line style, using points.");
				ostr << strPointStyle;
				break;
			case STYLE_POINTS:
				ostr << strPointStyle;
				break;
		}

		if(obj.oiColor)
		{
			unsigned int iColor = *obj.oiColor;
			char chFill = ostr.fill();
			ostr << std::setfill('0');
			ostr << "linecolor rgb \"#" << std::hex
				<< std::setw(2) << ((iColor & 0xff0000) >> 16)
				<< std::setw(2) << ((iColor & 0x00ff00) >> 8)
				<< std::setw(2) << ((iColor & 0x0000ff))
				<< "\" " << std::dec;
			ostr << std::setfill(chFill);
		}

		if(obj.strLegend != "")
		{
			m_bHasLegend = 1;
			ostr << " title \"" << obj.strLegend << "\" ";
		}
		else
		{
			ostr << " notitle ";
		}

		if(&obj != &(*m_vecObjs.rbegin()))
			ostr << ", ";
	}
	ostr << "\n";

	for(const PlotObj<t_real>& obj : m_vecObjs)
	{
		std::string strTab = BuildTable(obj.vecX, obj.vecY, obj.vecErrY, obj.vecErrX);
		ostr << strTab;
	}

	return ostr.str();
}


template<class t_real>
std::string GnuPlot<t_real>::BuildTable(const std::vector<t_real>& vecX, const std::vector<t_real>& vecY,
	const std::vector<t_real>& vecYErr, const std::vector<t_real>& vecXErr)
{
	std::ostringstream ostr;
	ostr.precision(m_iPrec);

	const std::size_t iSize = std::min(vecX.size(), vecY.size());
	const bool bHasXErr = (vecXErr.size() != 0);
	const bool bHasYErr = (vecYErr.size() != 0);

	for(std::size_t iDat=0; iDat<iSize; ++iDat)
	{
		ostr << vecX[iDat] << " " << vecY[iDat];

		if(bHasXErr) ostr << " " << vecXErr[iDat];
		if(bHasYErr) ostr << " " << vecYErr[iDat];

		ostr << "\n";
	}
	ostr << "end\n";

	return ostr.str();
}


template<class t_real>
void GnuPlot<t_real>::SetXLabel(const char* pcLab)
{
	if(!IsReady()) return;

	m_pProc->GetOstr() << "set xlabel \"" << pcLab << "\"\n";
}


template<class t_real>
void GnuPlot<t_real>::SetYLabel(const char* pcLab)
{
	if(!IsReady()) return;

	m_pProc->GetOstr() << "set ylabel \"" << pcLab << "\"\n";
}


template<class t_real>
void GnuPlot<t_real>::SetTitle(const char* pcTitle)
{
	if(!IsReady()) return;

	m_pProc->GetOstr() << "set title \"" << pcTitle << "\"\n";
}


template<class t_real>
void GnuPlot<t_real>::SetXRange(t_real dMin, t_real dMax)
{
	if(!IsReady()) return;

	m_pProc->GetOstr() << "set xrange [" << dMin << ":" << dMax << "]\n";
	m_pProc->GetOstr().flush();
}


template<class t_real>
void GnuPlot<t_real>::SetYRange(t_real dMin, t_real dMax)
{
	if(!IsReady()) return;

	m_pProc->GetOstr() << "set yrange [" << dMin << ":" << dMax << "]\n";
	m_pProc->GetOstr().flush();
}


template<class t_real>
void GnuPlot<t_real>::SetLogX(t_real tBase)
{
	if(!IsReady()) return;

	if(tBase >= 0.)
		m_pProc->GetOstr() << "set logscale x " << tBase << "\n";
	else
		m_pProc->GetOstr() << "unset logscale x\n";
	m_pProc->GetOstr().flush();
}


template<class t_real>
void GnuPlot<t_real>::SetLogY(t_real tBase)
{
	if(!IsReady()) return;

	if(tBase >= 0.)
		m_pProc->GetOstr() << "set logscale y " << tBase << "\n";
	else
		m_pProc->GetOstr() << "unset logscale y\n";
	m_pProc->GetOstr().flush();
}


template<class t_real>
void GnuPlot<t_real>::SetGrid(bool bOn)
{
	if(!IsReady()) return;

	if(bOn)
		m_pProc->GetOstr() << "set grid\n";
	else
		m_pProc->GetOstr() << "unset grid\n";

	m_pProc->GetOstr().flush();
}


template<class t_real>
void GnuPlot<t_real>::ClearArrows()
{
	if(!IsReady()) return;

	m_pProc->GetOstr() << "unset arrow\n";
	m_pProc->GetOstr().flush();
}


template<class t_real>
void GnuPlot<t_real>::AddArrow(t_real dX0, t_real dY0, t_real dX1, t_real dY1, bool bHead)
{
	if(!IsReady()) return;

	m_pProc->GetOstr() << "set arrow from " << dX0 << "," << dY0
			<< " to " << dX1 << "," << dY1;
	if(!bHead)
		m_pProc->GetOstr() << " nohead";

	m_pProc->GetOstr() << "\n";
	m_pProc->GetOstr().flush();
}


template<class t_real>
void GnuPlot<t_real>::SetColorBarRange(t_real dMin, t_real dMax, bool bCyclic)
{
	if(!IsReady()) return;

	m_pProc->GetOstr() << "set cbrange [" << dMin << ":" << dMax << "]\n";

	if(bCyclic)
		m_pProc->GetOstr() << "set palette defined (0 \"#0000ff\", 0.33333 \"#ff0000\", 0.66666 \"#ff9900\", 1 \"#0000ff\")\n";
	else
		m_pProc->GetOstr() << "set palette defined (0 \"#0000ff\", 1 \"#ff0000\")\n";

	m_pProc->GetOstr().flush();
}


template<class t_real>
bool GnuPlot<t_real>::IsReady() const { return m_pProc && m_pProc->IsReady(); }

template<class t_real>
std::ostream& GnuPlot<t_real>::GetStream() { return m_pProc->GetOstr(); }

}

#endif
