/**
 * Loads instrument-specific data files
 * @author Tobias Weber <tobias.weber@tum.de>
 * @date feb-2015
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

#ifndef __LOADINSTR_H__
#define __LOADINSTR_H__

#include <unordered_map>
#include <map>
#include <vector>
#include <array>
#include <iostream>
#include "../string/string.h"
#include "loaddat.h"


namespace tl{

// interface for instrument-specific data files
template<class _t_real = double>
class FileInstrBase
{
	public:
		using t_real = _t_real;

		typedef std::unordered_map<std::string, std::string> t_mapParams;
		typedef std::vector<std::string> t_vecColNames;
		typedef std::vector<t_real> t_vecVals;
		typedef std::vector<t_vecVals> t_vecDat;

	protected:
		void RenameDuplicateCols();

		std::array<t_real, 5> GetScanHKLKiKf(const char* pcH, const char* pcK,
			const char* pcL, const char* pcE, std::size_t i) const;

	public:
		FileInstrBase() = default;
		virtual ~FileInstrBase() = default;

		virtual bool Load(const char* pcFile) = 0;

		virtual std::array<t_real, 3> GetSampleLattice() const = 0;
		virtual std::array<t_real, 3> GetSampleAngles() const = 0;
		virtual std::array<t_real, 2> GetMonoAnaD() const = 0;

		virtual std::array<bool, 3> GetScatterSenses() const = 0;
		virtual std::array<t_real, 3> GetScatterPlane0() const = 0;
		virtual std::array<t_real, 3> GetScatterPlane1() const = 0;

		virtual t_real GetKFix() const = 0;
		virtual bool IsKiFixed() const = 0;

		virtual const t_vecVals& GetCol(const std::string& strName, std::size_t *pIdx=0) const = 0;
		virtual t_vecVals& GetCol(const std::string& strName, std::size_t *pIdx=0) = 0;

		virtual std::array<t_real, 4> GetPosHKLE() const = 0;	// zero pos.

		virtual std::size_t GetScanCount() const = 0;
		virtual std::array<t_real, 5> GetScanHKLKiKf(std::size_t i) const = 0;
		virtual bool MergeWith(const FileInstrBase<t_real>* pDat);
		virtual void SmoothData(const std::string& strCol, t_real dEps, bool bIterate=1);

		virtual std::string GetTitle() const = 0;
		virtual std::string GetUser() const = 0;
		virtual std::string GetLocalContact() const = 0;
		virtual std::string GetScanNumber() const = 0;
		virtual std::string GetSampleName() const = 0;
		virtual std::string GetSpacegroup() const = 0;
		virtual std::string GetTimestamp() const = 0;

		virtual const t_vecDat& GetData() const = 0;
		virtual t_vecDat& GetData() = 0;
		virtual const t_vecColNames& GetColNames() const = 0;
		virtual const t_mapParams& GetAllParams() const = 0;

		virtual std::vector<std::string> GetScannedVars() const = 0;
		virtual std::string GetCountVar() const = 0;
		virtual std::string GetMonVar() const = 0;

		virtual std::string GetScanCommand() const = 0;

		// polarisation stuff
		virtual void ParsePolData();
		virtual std::size_t NumPolChannels() const;
		virtual const std::vector<std::array<t_real, 6>>& GetPolStates() const;
		virtual void SetPolNames(const char* pVec1, const char* pVec2,
			const char* pCur1, const char* pCur2);	// spherical PA
		virtual void SetLinPolNames(const char* pFlip1, const char* pFlip2,
			const char* pXYZ);	// linear PA

	public:
		virtual bool MatchColumn(const std::string& strRegex,
			std::string& strColName, bool bSortByCounts=0, bool bFilterEmpty=1) const;

		static FileInstrBase<t_real>* LoadInstr(const char* pcFile);
};


// psi & ill files
template<class _t_real = double>
class FilePsi : public FileInstrBase<_t_real>
{
	public:
		using t_real = _t_real;
		using t_mapParams = typename FileInstrBase<t_real>::t_mapParams;
		using t_vecColNames = typename FileInstrBase<t_real>::t_vecColNames;
		using t_vecVals = typename FileInstrBase<t_real>::t_vecVals;
		using t_vecDat = typename FileInstrBase<t_real>::t_vecDat;

		// internal parameters in m_mapParams
		typedef std::map<std::string, t_real> t_mapIParams;

	protected:
		t_mapParams m_mapParams;
		t_mapIParams m_mapParameters, m_mapZeros, m_mapVariables, m_mapPosHkl, m_mapScanSteps;
		t_vecColNames m_vecColNames;
		t_vecDat m_vecData;

		// automatically look and parse for polarisation data
		bool m_bAutoParsePol = false;

		// incoming and outgoing polarisation states
		std::vector<std::array<t_real, 6>> m_vecPolStates;

		// instrument-specific device names
		std::string m_strPolVec1 {"p1"}, m_strPolVec2 {"p2"};
		std::string m_strPolCur1 {"i1"}, m_strPolCur2 {"i2"};

		std::string m_strXYZ {"hx"};
		std::string m_strFlip1 {"f1"}, m_strFlip2 {"f2"};


	protected:
		std::string ReadData(std::istream& istr);
		std::string ReadMultiData(std::istream& istr);
		void GetInternalParams(const std::string& strAll, t_mapIParams& mapPara, bool fix_broken = false);

	public:
		FilePsi() = default;
		virtual ~FilePsi() = default;

		virtual bool Load(const char* pcFile) override;

		void PrintParams(std::ostream& ostr) const;
		const t_mapParams& GetParams() const { return m_mapParams; }

		const std::string& GetColName(std::size_t iCol) const { return m_vecColNames[iCol]; }
		std::size_t GetColCount() const { return m_vecColNames.size(); }

		bool HasCol(const std::string& strName) const;
		const t_vecVals& GetCol(std::size_t iCol) const { return m_vecData[iCol]; }
		virtual const t_vecVals& GetCol(const std::string& strName, std::size_t *pIdx=0) const override;
		virtual t_vecVals& GetCol(const std::string& strName, std::size_t *pIdx=0) override;

	public:
		virtual std::array<t_real, 3> GetSampleLattice() const override;
		virtual std::array<t_real, 3> GetSampleAngles() const override;
		virtual std::array<t_real, 2> GetMonoAnaD() const override;

		virtual std::array<bool, 3> GetScatterSenses() const override;
		virtual std::array<t_real, 3> GetScatterPlane0() const override;
		virtual std::array<t_real, 3> GetScatterPlane1() const override;

		virtual t_real GetKFix() const override;
		virtual bool IsKiFixed() const override;

		virtual std::array<t_real, 4> GetPosHKLE() const override;	// zero pos.
		std::array<t_real, 4> GetDeltaHKLE() const;	// scan steps

		virtual std::size_t GetScanCount() const override;
		virtual std::array<t_real, 5> GetScanHKLKiKf(std::size_t i) const override;
		virtual bool MergeWith(const FileInstrBase<t_real>* pDat) override;

		virtual std::string GetTitle() const override;
		virtual std::string GetUser() const override;
		virtual std::string GetLocalContact() const override;
		virtual std::string GetScanNumber() const override;
		virtual std::string GetSampleName() const override;
		virtual std::string GetSpacegroup() const override;
		virtual std::string GetTimestamp() const override;

		virtual const t_vecDat& GetData() const override { return m_vecData; }
		virtual t_vecDat& GetData() override { return m_vecData; }
		virtual const t_vecColNames& GetColNames() const override { return m_vecColNames; }
		virtual const t_mapParams& GetAllParams() const override { return m_mapParams; }

		virtual std::vector<std::string> GetScannedVars() const override;
		virtual std::string GetCountVar() const override;
		virtual std::string GetMonVar() const override;

		virtual std::string GetScanCommand() const override;

		virtual std::size_t NumPolChannels() const override;

	public:
		void SetAutoParsePolData(bool b) { m_bAutoParsePol = b; }
		virtual void ParsePolData() override;

		virtual const std::vector<std::array<t_real, 6>>& GetPolStates() const override
		{
			return m_vecPolStates;
		}

		virtual void SetPolNames(const char* pVec1, const char* pVec2,
			const char* pCur1, const char* pCur2) override
		{
			m_strPolVec1 = pVec1; m_strPolVec2 = pVec2;
			m_strPolCur1 = pCur1; m_strPolCur2 = pCur2;
		}

		virtual void SetLinPolNames(const char* pFlip1, const char* pFlip2,
			const char* pXYZ) override
		{
			m_strFlip1 = pFlip1; m_strFlip2 = pFlip2;
			m_strXYZ = pXYZ;
		}
};


// frm/nicos files
template<class _t_real = double>
class FileFrm : public FileInstrBase<_t_real>
{
	public:
		using t_real = _t_real;
		using t_mapParams = typename FileInstrBase<t_real>::t_mapParams;
		using t_vecColNames = typename FileInstrBase<t_real>::t_vecColNames;
		using t_vecVals = typename FileInstrBase<t_real>::t_vecVals;
		using t_vecDat = typename FileInstrBase<t_real>::t_vecDat;

	protected:
		t_mapParams m_mapParams;
		t_vecColNames m_vecQuantities, m_vecUnits;
		t_vecDat m_vecData;
		std::string m_strInstrIdent;

	public:
		FileFrm() = default;
		virtual ~FileFrm() = default;

	protected:
		void ReadHeader(std::istream& istr);
		void ReadData(std::istream& istr);

	public:
		virtual bool Load(const char* pcFile) override;

		virtual std::array<t_real, 3> GetSampleLattice() const override;
		virtual std::array<t_real, 3> GetSampleAngles() const override;
		virtual std::array<t_real, 2> GetMonoAnaD() const override;

		virtual std::array<bool, 3> GetScatterSenses() const override;
		virtual std::array<t_real, 3> GetScatterPlane0() const override;
		virtual std::array<t_real, 3> GetScatterPlane1() const override;

		virtual std::array<t_real, 4> GetPosHKLE() const override;	// zero pos.

		virtual t_real GetKFix() const override;
		virtual bool IsKiFixed() const override;

		virtual std::size_t GetScanCount() const override;
		virtual std::array<t_real, 5> GetScanHKLKiKf(std::size_t i) const override;
		virtual bool MergeWith(const FileInstrBase<t_real>* pDat) override;

		virtual const t_vecVals& GetCol(const std::string& strName, std::size_t *pIdx=0) const override;
		virtual t_vecVals& GetCol(const std::string& strName, std::size_t *pIdx=0) override;

		virtual std::string GetTitle() const override;
		virtual std::string GetUser() const override;
		virtual std::string GetLocalContact() const override;
		virtual std::string GetScanNumber() const override;
		virtual std::string GetSampleName() const override;
		virtual std::string GetSpacegroup() const override;
		virtual std::string GetTimestamp() const override;

		virtual const t_vecDat& GetData() const override { return m_vecData; }
		virtual t_vecDat& GetData() override { return m_vecData; }
		virtual const t_vecColNames& GetColNames() const override { return m_vecQuantities; }
		virtual const t_mapParams& GetAllParams() const override { return m_mapParams; }

		virtual std::vector<std::string> GetScannedVars() const override;
		virtual std::string GetCountVar() const override;
		virtual std::string GetMonVar() const override;

		virtual std::string GetScanCommand() const override;
};


// macs files
template<class _t_real = double>
class FileMacs : public FileInstrBase<_t_real>
{
	public:
		using t_real = _t_real;
		using t_mapParams = typename FileInstrBase<t_real>::t_mapParams;
		using t_vecColNames = typename FileInstrBase<t_real>::t_vecColNames;
		using t_vecVals = typename FileInstrBase<t_real>::t_vecVals;
		using t_vecDat = typename FileInstrBase<t_real>::t_vecDat;

	protected:
		t_mapParams m_mapParams;
		t_vecColNames m_vecQuantities;
		t_vecDat m_vecData;

	public:
		FileMacs() = default;
		virtual ~FileMacs() = default;

	protected:
		void ReadHeader(std::istream& istr);
		void ReadData(std::istream& istr);

	public:
		virtual bool Load(const char* pcFile) override;

		virtual std::array<t_real, 3> GetSampleLattice() const override;
		virtual std::array<t_real, 3> GetSampleAngles() const override;
		virtual std::array<t_real, 2> GetMonoAnaD() const override;

		virtual std::array<bool, 3> GetScatterSenses() const override;
		virtual std::array<t_real, 3> GetScatterPlane0() const override;
		virtual std::array<t_real, 3> GetScatterPlane1() const override;

		virtual std::array<t_real, 4> GetPosHKLE() const override;	// zero pos.

		virtual t_real GetKFix() const override;
		virtual bool IsKiFixed() const override;

		virtual std::size_t GetScanCount() const override;
		virtual std::array<t_real, 5> GetScanHKLKiKf(std::size_t i) const override;
		virtual bool MergeWith(const FileInstrBase<t_real>* pDat) override;

		virtual const t_vecVals& GetCol(const std::string& strName, std::size_t *pIdx=0) const override;
		virtual t_vecVals& GetCol(const std::string& strName, std::size_t *pIdx=0) override;

		virtual std::string GetTitle() const override;
		virtual std::string GetUser() const override;
		virtual std::string GetLocalContact() const override;
		virtual std::string GetScanNumber() const override;
		virtual std::string GetSampleName() const override;
		virtual std::string GetSpacegroup() const override;
		virtual std::string GetTimestamp() const override;

		virtual const t_vecDat& GetData() const override { return m_vecData; }
		virtual t_vecDat& GetData() override { return m_vecData; }
		virtual const t_vecColNames& GetColNames() const override { return m_vecQuantities; }
		virtual const t_mapParams& GetAllParams() const override { return m_mapParams; }

		virtual std::vector<std::string> GetScannedVars() const override;
		virtual std::string GetCountVar() const override;
		virtual std::string GetMonVar() const override;

		virtual std::string GetScanCommand() const override;
};


// trisp files
template<class _t_real = double>
class FileTrisp : public FileInstrBase<_t_real>
{
	public:
		using t_real = _t_real;
		using t_mapParams = typename FileInstrBase<t_real>::t_mapParams;
		using t_vecColNames = typename FileInstrBase<t_real>::t_vecColNames;
		using t_vecVals = typename FileInstrBase<t_real>::t_vecVals;
		using t_vecDat = typename FileInstrBase<t_real>::t_vecDat;

	protected:
		t_mapParams m_mapParams;
		t_vecColNames m_vecQuantities;
		t_vecDat m_vecData;

	public:
		FileTrisp() = default;
		virtual ~FileTrisp() = default;

	protected:
		void ReadHeader(std::istream& istr);
		void ReadData(std::istream& istr);

	public:
		virtual bool Load(const char* pcFile) override;

		virtual std::array<t_real, 3> GetSampleLattice() const override;
		virtual std::array<t_real, 3> GetSampleAngles() const override;
		virtual std::array<t_real, 2> GetMonoAnaD() const override;

		virtual std::array<bool, 3> GetScatterSenses() const override;
		virtual std::array<t_real, 3> GetScatterPlane0() const override;
		virtual std::array<t_real, 3> GetScatterPlane1() const override;

		virtual std::array<t_real, 4> GetPosHKLE() const override;	// zero pos.

		virtual t_real GetKFix() const override;
		virtual bool IsKiFixed() const override;

		virtual std::size_t GetScanCount() const override;
		virtual std::array<t_real, 5> GetScanHKLKiKf(std::size_t i) const override;
		virtual bool MergeWith(const FileInstrBase<t_real>* pDat) override;

		virtual const t_vecVals& GetCol(const std::string& strName, std::size_t *pIdx=0) const override;
		virtual t_vecVals& GetCol(const std::string& strName, std::size_t *pIdx=0) override;

		virtual std::string GetTitle() const override;
		virtual std::string GetUser() const override;
		virtual std::string GetLocalContact() const override;
		virtual std::string GetScanNumber() const override;
		virtual std::string GetSampleName() const override;
		virtual std::string GetSpacegroup() const override;
		virtual std::string GetTimestamp() const override;

		virtual const t_vecDat& GetData() const override { return m_vecData; }
		virtual t_vecDat& GetData() override { return m_vecData; }
		virtual const t_vecColNames& GetColNames() const override { return m_vecQuantities; }
		virtual const t_mapParams& GetAllParams() const override { return m_mapParams; }

		virtual std::vector<std::string> GetScannedVars() const override;
		virtual std::string GetCountVar() const override;
		virtual std::string GetMonVar() const override;

		virtual std::string GetScanCommand() const override;
};


// raw data files
template<class _t_real = double>
class FileRaw : public FileInstrBase<_t_real>
{
	public:
		using t_real = _t_real;
		using t_mapParams = typename FileInstrBase<t_real>::t_mapParams;
		using t_vecColNames = typename FileInstrBase<t_real>::t_vecColNames;
		using t_vecVals = typename FileInstrBase<t_real>::t_vecVals;
		using t_vecDat = typename FileInstrBase<t_real>::t_vecDat;

	protected:
		DatFile<t_real, char> m_dat;
		t_vecColNames m_vecCols;

	public:
		FileRaw() = default;
		virtual ~FileRaw() = default;

	public:
		virtual bool Load(const char* pcFile) override;

		virtual std::array<t_real, 3> GetSampleLattice() const override;
		virtual std::array<t_real, 3> GetSampleAngles() const override;
		virtual std::array<t_real, 2> GetMonoAnaD() const override;

		virtual std::array<bool, 3> GetScatterSenses() const override;
		virtual std::array<t_real, 3> GetScatterPlane0() const override;
		virtual std::array<t_real, 3> GetScatterPlane1() const override;

		virtual std::array<t_real, 4> GetPosHKLE() const override;	// zero pos.

		virtual t_real GetKFix() const override;
		virtual bool IsKiFixed() const override;

		virtual std::size_t GetScanCount() const override;
		virtual std::array<t_real, 5> GetScanHKLKiKf(std::size_t i) const override;
		virtual bool MergeWith(const FileInstrBase<t_real>* pDat) override;

		virtual const t_vecVals& GetCol(const std::string& strName, std::size_t *pIdx=0) const override;
		virtual t_vecVals& GetCol(const std::string& strName, std::size_t *pIdx=0) override;

		virtual std::string GetTitle() const override;
		virtual std::string GetUser() const override;
		virtual std::string GetLocalContact() const override;
		virtual std::string GetScanNumber() const override;
		virtual std::string GetSampleName() const override;
		virtual std::string GetSpacegroup() const override;
		virtual std::string GetTimestamp() const override;

		virtual const t_vecDat& GetData() const override;
		virtual t_vecDat& GetData() override;
		virtual const t_vecColNames& GetColNames() const override;
		virtual const t_mapParams& GetAllParams() const override;

		virtual std::vector<std::string> GetScannedVars() const override;
		virtual std::string GetCountVar() const override;
		virtual std::string GetMonVar() const override;

		virtual std::string GetScanCommand() const override;
};


#ifdef USE_HDF5
// hdf5 data files
template<class _t_real = double>
class FileH5 : public FileInstrBase<_t_real>
{
public:
	using t_real = _t_real;
	using t_mapParams = typename FileInstrBase<t_real>::t_mapParams;
	using t_vecColNames = typename FileInstrBase<t_real>::t_vecColNames;
	using t_vecVals = typename FileInstrBase<t_real>::t_vecVals;
	using t_vecDat = typename FileInstrBase<t_real>::t_vecDat;

protected:
	t_vecDat m_data;
	t_vecColNames m_vecCols;
	t_mapParams m_params;

	std::string m_title, m_username, m_localname, m_timestamp;
	int m_scannumber = 0;
	std::string m_scancommand;

	std::array<t_real, 3> m_lattice, m_angles, m_plane[2];
	std::array<t_real, 4> m_initialpos;

	std::array<t_real, 2> m_dspacings;
	std::array<bool, 3> m_senses;

	t_real m_kfix = 0.;
	bool m_iskifixed = false;

	std::vector<std::string> m_scanned_vars;

public:
	FileH5() = default;
	virtual ~FileH5() = default;

public:
	virtual bool Load(const char* pcFile) override;

	virtual std::array<t_real, 3> GetSampleLattice() const override;
	virtual std::array<t_real, 3> GetSampleAngles() const override;
	virtual std::array<t_real, 2> GetMonoAnaD() const override;

	virtual std::array<bool, 3> GetScatterSenses() const override;
	virtual std::array<t_real, 3> GetScatterPlane0() const override;
	virtual std::array<t_real, 3> GetScatterPlane1() const override;

	virtual std::array<t_real, 4> GetPosHKLE() const override;

	virtual t_real GetKFix() const override;
	virtual bool IsKiFixed() const override;

	virtual std::size_t GetScanCount() const override;
	virtual std::array<t_real, 5> GetScanHKLKiKf(std::size_t i) const override;
	virtual bool MergeWith(const FileInstrBase<t_real>* pDat) override;

	std::size_t GetColCount() const { return m_vecCols.size(); }
	virtual const t_vecVals& GetCol(const std::string& strName, std::size_t *pIdx=0) const override;
	virtual t_vecVals& GetCol(const std::string& strName, std::size_t *pIdx=0) override;

	virtual std::string GetTitle() const override;
	virtual std::string GetUser() const override;
	virtual std::string GetLocalContact() const override;
	virtual std::string GetScanNumber() const override;
	virtual std::string GetSampleName() const override;
	virtual std::string GetSpacegroup() const override;
	virtual std::string GetTimestamp() const override;

	virtual const t_vecDat& GetData() const override;
	virtual t_vecDat& GetData() override;
	virtual const t_vecColNames& GetColNames() const override;
	virtual const t_mapParams& GetAllParams() const override;

	virtual std::vector<std::string> GetScannedVars() const override;
	virtual std::string GetCountVar() const override;
	virtual std::string GetMonVar() const override;

	virtual std::string GetScanCommand() const override;
};
#endif

}

#ifdef TLIBS_INC_HDR_IMPLS
	#include "loadinstr_impl.h"
#endif

#endif
