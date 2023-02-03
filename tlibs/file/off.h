/**
 * off file handling
 * @author Tobias Weber <tobias.weber@tum.de>
 * @date 23-sep-17
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

#ifndef __OFF_FILES_H__
#define __OFF_FILES_H__

#include <ostream>
#include <fstream>
#include <vector>
#include <list>
#include <unordered_set>
#include <string>

#include "../math/linalg.h"
#include "../string/string.h"


namespace tl{

enum class Off3dReadOp
{
	READ_MAGIC,
	READ_HEADER,
	READ_VERTICES,
	READ_POLYS,
	READ_FINISHED
};


template<class t_real = double>
class Off3d
{
public:
	template<class... Args> using t_vec_bare = ublas::vector<Args...>;
	using t_vec = t_vec_bare<t_real>;

protected:
	std::list<t_vec> m_lstVertices;
	std::list<std::list<std::size_t>> m_lstPolys;
	std::string m_strComment;

public:
	Off3d() = default;
	~Off3d() = default;

	void SetComment(const std::string& str) { m_strComment = str; }

	// -------------------------------------------------------------------------

	/**
	 * remove vertex with index iIdx and replace all references with iIdxReplace
	 */
	void RemoveVertex(std::size_t iIdx, std::size_t iIdxReplace)
	{
		typename decltype(m_lstVertices)::iterator iter = m_lstVertices.begin();
		std::advance(iter, iIdx);
		m_lstVertices.erase(iter);

		for(std::list<std::size_t>& lstPoly : m_lstPolys)
		{
			for(std::size_t& iIdxPoly : lstPoly)
			{
				// replace removed index with given new one
				if(iIdxPoly == iIdx)
					iIdxPoly = iIdxReplace;
				// decrease following indices
				else if(iIdxPoly > iIdx)
					--iIdxPoly;
			}
		}
	}


	void Optimise(t_real eps = std::numeric_limits<t_real>::epsilon())
	{
		// nothing to do
		if(m_lstVertices.size() == 0)
			return;


		// -------------------------------------
		// remove unreferenced vertices
		std::unordered_set<std::size_t> setIdx;
		for(std::list<std::size_t>& lstPoly : m_lstPolys)
		{
			for(std::size_t iIdxPoly : lstPoly)
				setIdx.insert(iIdxPoly);
		}
		for(std::ptrdiff_t iVert = m_lstVertices.size()-1; iVert >= 0; --iVert)
		{
			if(setIdx.find(iVert) == setIdx.end())
				RemoveVertex(std::size_t(iVert), 0);
		}
		// -------------------------------------


		// -------------------------------------
		// remove duplicate vertices (within epsilon range)
		std::vector<std::pair<std::size_t, std::size_t>> vecVertsToRemove;

		for(std::size_t iPt0 = 0; iPt0 < m_lstVertices.size(); ++iPt0)
		{
			const t_vec& vec0 = *std::next(m_lstVertices.begin(), iPt0);

			for(std::size_t iPt1 = iPt0+1; iPt1 < m_lstVertices.size(); ++iPt1)
			{
				const t_vec& vec1 = *std::next(m_lstVertices.begin(), iPt1);

				if(vec_equal(vec0, vec1, eps))
					vecVertsToRemove.push_back(std::make_pair(iPt1, iPt0));
			}
		}

		for(const std::pair<std::size_t, std::size_t>& pair : vecVertsToRemove)
			RemoveVertex(pair.first, pair.second);
		// -------------------------------------
	}

	// -------------------------------------------------------------------------

	/**
	 * load file from stream
	 */
	bool Read(std::istream& istr)
	{
		Off3dReadOp op = Off3dReadOp::READ_MAGIC;
		std::size_t iNumVerts = 0, iNumPolys = 0;
		std::size_t iCurVert = 0, iCurPoly = 0;

		while(!istr.eof())
		{
			std::string strLine;
			std::getline(istr, strLine);
			if(strLine == "")
				continue;
			if(strLine[0] == '#')
				continue;

			// all vertices read? -> read polys next
			if(op == Off3dReadOp::READ_VERTICES && iCurVert >= iNumVerts)
				op = Off3dReadOp::READ_POLYS;
			// all polygons read? -> finish
			if(op == Off3dReadOp::READ_POLYS && iCurPoly >= iNumPolys)
				op = Off3dReadOp::READ_FINISHED;


			if(op == Off3dReadOp::READ_MAGIC)
			{
				if(strLine == "OFF")
					op = Off3dReadOp::READ_HEADER;
			}
			else if(op == Off3dReadOp::READ_HEADER)
			{
				std::istringstream istrstr(strLine);
				istrstr >> iNumVerts >> iNumPolys;
				op = Off3dReadOp::READ_VERTICES;
			}
			else if(op == Off3dReadOp::READ_VERTICES)
			{
				std::vector<t_real> vec;
				get_tokens<t_real, std::string, decltype(vec)>(strLine, " \t", vec);
				t_vec vert = convert_vec_full<t_real, t_real, std::vector, t_vec_bare>(vec);

				m_lstVertices.emplace_back(std::move(vert));
				++iCurVert;
			}
			else if(op == Off3dReadOp::READ_POLYS)
			{
				std::list<std::size_t> lstIdx;
				get_tokens<std::size_t, std::string, decltype(lstIdx)>(strLine, " \t", lstIdx);
				if(lstIdx.size() == 0)
					continue;

				std::size_t iLen = lstIdx.front();
				lstIdx.pop_front();
				assert(lstIdx.size() == iLen);

				m_lstPolys.push_back(lstIdx);
				++iCurPoly;
			}
			else if(op == Off3dReadOp::READ_FINISHED)
			{
				break;
			}
		}

		if(iCurPoly == iNumPolys && iCurVert == iNumVerts)
			op = Off3dReadOp::READ_FINISHED;

		return (op == Off3dReadOp::READ_FINISHED);
	}


	/**
	 * write file to stream
	 */
	bool Write(std::ostream& ostr) const
	{
		// Magic & header
		ostr << "OFF\n\n";

		if(m_strComment != "")
			ostr << "# " << m_strComment << "\n\n";

		ostr << m_lstVertices.size() << " "
			<< m_lstPolys.size() << " "
			<< "0\n\n";


		// write vertices
		ostr << "# vertices\n";

		for(const t_vec& vert : m_lstVertices)
		{
			for(std::size_t i = 0; i < vert.size(); ++i)
			{
				ostr << vert[i];
				if(i != vert.size()-1)
					ostr << " ";
			}
			ostr << "\n";
		}
		ostr << "\n";


		// write polygons
		ostr << "# polygons\n";

		for(const std::list<std::size_t>& lstPoly : m_lstPolys)
		{
			ostr << lstPoly.size() << " ";
	
			for(std::list<std::size_t>::const_iterator iter = lstPoly.begin();
				iter != lstPoly.end(); ++iter)
			{
				ostr << *iter;
				if(std::next(iter) != lstPoly.end())
					ostr << " ";
			}
			ostr << "\n";
		}
		ostr << "\n";


		ostr.flush();
		return true;
	}


	/**
	 * save to file
	 */
	bool Save(const char* pcFile) const
	{
		std::ofstream ofstr(pcFile);
		if(!ofstr)
			return false;

		return Write(ofstr);
	}


	/**
	 * load from file
	 */
	bool Load(const char* pcFile)
	{
		std::ifstream ifstr(pcFile);
		if(!ifstr)
			return false;

		return Read(ifstr);
	}
};

}

#endif
