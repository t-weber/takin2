/**
 * simplified x3d file handling
 * @author Tobias Weber <tobias.weber@tum.de>
 * @date nov-2015
 * @license GPLv2 or GPLv3
 */

#ifndef __X3D_FILES_H__
#define __X3D_FILES_H__

#include <ostream>
#include <fstream>
#include <vector>
#include <string>

#include "../math/linalg.h"
#include "../math/quat.h"

namespace tl{


template<class t_real = double> class X3dElem;

template<class t_real>
class X3dElem
{
	public:
		using t_vec = ublas::vector<t_real>;
		using t_mat = ublas::matrix<t_real>;
		using t_quat = math::quaternion<t_real>;

	protected:
		static void WriteColor(std::ostream& ostr, const t_vec& vecColor)
		{
			if(vecColor.size() >= 3)
			{
				ostr << "<Appearance>\n";
					ostr << "<Material diffuseColor=\""
						<< vecColor[0] << " " << vecColor[1] << " " << vecColor[2]
						<< "\" ";
				if(vecColor.size() >= 4)
					ostr << "transparency=\"" << vecColor[3] << "\" ";
				ostr << "/>\n";
				ostr << "</Appearance>\n";
			}

		}

	protected:
		std::vector<X3dElem<t_real>*> m_vecChildren;

	public:
		X3dElem() = default;
		virtual ~X3dElem()
		{
			for(X3dElem* pElem : m_vecChildren)
				delete pElem;
			m_vecChildren.clear();
		}

		virtual void Write(std::ostream& ostr) const
		{
			for(const X3dElem *pElem : m_vecChildren)
				pElem->Write(ostr);
		}

		virtual void AddChild(X3dElem* pElem)
		{
			if(pElem)
				m_vecChildren.push_back(pElem);
		}
};


// -----------------------------------------------------------------------------


template<class t_real = double>
class X3dScene : public X3dElem<t_real>
{
	public:
		X3dScene() = default;
		virtual ~X3dScene() = default;

		virtual void Write(std::ostream& ostr) const
		{
			ostr << "<Scene>\n";
			X3dElem<t_real>::Write(ostr);
			ostr << "</Scene>\n";
		}
};


// -----------------------------------------------------------------------------


template<class t_real = double>
class X3dTrafo : public X3dElem<t_real>
{
	public:
		using typename X3dElem<t_real>::t_vec;
		using typename X3dElem<t_real>::t_mat;
		using typename X3dElem<t_real>::t_quat;

	protected:
		t_vec m_vecTrans;
		t_vec m_vecScale;
		t_quat m_quatRot;
		bool m_bHasRot = 0;

	public:
		X3dTrafo() = default;
		virtual ~X3dTrafo() = default;

		virtual void Write(std::ostream& ostr) const override
		{
			ostr << "<Transform ";

			if(m_vecTrans.size() >= 3)
			{
				ostr << "translation=\""
					<< m_vecTrans[0] << " " << m_vecTrans[1] << " " << m_vecTrans[2]
					<< "\" ";
			}
			if(m_vecScale.size() >= 3)
			{
				ostr << "scale=\""
					<< m_vecScale[0] << " " << m_vecScale[1] << " " << m_vecScale[2]
					<< "\" ";
			}
			if(m_bHasRot)
			{
				ostr << "rotation=\""
					<< m_quatRot.R_component_1() << " " << m_quatRot.R_component_2() << " "
					<< m_quatRot.R_component_3() << " " << m_quatRot.R_component_4()
					<< "\" ";
			}

			ostr << ">\n";

			X3dElem<t_real>::Write(ostr);

			ostr << "</Transform>\n";
		}

		void SetTrans(const t_vec& vec) { m_vecTrans = vec; }
		void SetScale(const t_vec& vec) { m_vecScale = vec; }
		void SetRot(const t_quat& quat) { m_quatRot = quat; m_bHasRot = 1; }
};


// -----------------------------------------------------------------------------


// sphere
template<class t_real = double>
class X3dSphere : public X3dElem<t_real>
{
	public:
		using typename X3dElem<t_real>::t_vec;
		using typename X3dElem<t_real>::t_mat;
		using typename X3dElem<t_real>::t_quat;

	protected:
		t_real m_dRadius = 1.;
		t_vec m_vecColor;

	public:
		X3dSphere() = default;
		X3dSphere(t_real dRad) : m_dRadius(dRad) {}
		virtual ~X3dSphere() = default;

		virtual void Write(std::ostream& ostr) const override
		{
			ostr << "<Shape>\n";
				ostr << "<Sphere radius=\"" << m_dRadius << "\" />\n";

			X3dElem<t_real>::WriteColor(ostr, m_vecColor);
			ostr << "</Shape>\n";

			// TODO: child elements?
			//X3dElem<t_real>::Write(ostr);
		}

		void SetRadius(t_real dRad) { m_dRadius = dRad; }
		void SetColor(const t_vec& vecCol) { m_vecColor = vecCol; }
};

// -----------------------------------------------------------------------------

// cuboid
template<class t_real = double>
class X3dCube : public X3dElem<t_real>
{
	public:
		using typename X3dElem<t_real>::t_vec;
		using typename X3dElem<t_real>::t_mat;
		using typename X3dElem<t_real>::t_quat;

	protected:
		t_vec m_vecLength;
		t_vec m_vecColor;

	public:
		X3dCube() = default;
		X3dCube(t_real dW, t_real dH, t_real dL)
			: m_vecLength(tl::make_vec({dW, dH, dL})) {}
		virtual ~X3dCube() = default;

		virtual void Write(std::ostream& ostr) const override
		{
			ostr << "<Shape>\n";
				ostr << "<Box size=\""
					<< m_vecLength[0] << " "
					<< m_vecLength[1] << " "
					<< m_vecLength[2] << "\" />\n";

			X3dElem<t_real>::WriteColor(ostr, m_vecColor);
			ostr << "</Shape>\n";

			// TODO: child elements?
			//X3dElem<t_real>::Write(ostr);
		}

		void SetLengths(const t_vec& vec) { m_vecLength = vec; }
		void SetColor(const t_vec& vecCol) { m_vecColor = vecCol; }
};

// -----------------------------------------------------------------------------

// cylinder
template<class t_real = double>
class X3dCylinder : public X3dElem<t_real>
{
	public:
		using typename X3dElem<t_real>::t_vec;
		using typename X3dElem<t_real>::t_mat;
		using typename X3dElem<t_real>::t_quat;

	protected:
		t_real m_dRadius;
		t_real m_dHeight;
		t_vec m_vecColor;

	public:
		X3dCylinder() = default;
		X3dCylinder(t_real dRadius, t_real dHeight)
			: m_dRadius(dRadius), m_dHeight(dHeight) {}
		virtual ~X3dCylinder() = default;

		virtual void Write(std::ostream& ostr) const override
		{
			ostr << "<Shape>\n";
				ostr << "<Cylinder height=\"" << m_dHeight << "\""
					<< " radius=\"" << m_dRadius << "\" />\n";

			X3dElem<t_real>::WriteColor(ostr, m_vecColor);
			ostr << "</Shape>\n";

			// TODO: child elements?
			//X3dElem<t_real>::Write(ostr);
		}


		void SetRadius(const t_real dRad) { m_dRadius = dRad; }
		void SetHeight(const t_real dHeight) { m_dHeight = dHeight; }
		void SetColor(const t_vec& vecCol) { m_vecColor = vecCol; }
};

// -----------------------------------------------------------------------------

// polygon
template<class t_real = double>
class X3dPolygon : public X3dElem<t_real>
{
	public:
		using typename X3dElem<t_real>::t_vec;
		using typename X3dElem<t_real>::t_mat;
		using typename X3dElem<t_real>::t_quat;

	protected:
		std::vector<t_vec> m_vertices;
		t_vec m_vecColor;

	public:
		X3dPolygon() = default;
		virtual ~X3dPolygon() = default;

		virtual void Write(std::ostream& ostr) const override
		{
			std::ostringstream ostrVertices;
			ostrVertices << "<Coordinate point=\"";
			for(std::size_t iVert=0; iVert<m_vertices.size(); ++iVert)
			{
				const t_vec& vec = m_vertices[iVert];

				for(t_real d : vec)
					ostrVertices << d << " ";

				if(iVert+1 < m_vertices.size())
					ostrVertices << ", ";
			}
			ostrVertices << "\" />\n";


			ostr << "<Shape>\n";
				ostr << "<IndexedFaceSet ";
				ostr << "coordIndex=\"";
				for(std::size_t iVert=0; iVert<m_vertices.size(); ++iVert)
					ostr << iVert << ", ";
				ostr << "-1\"";
				ostr << ">\n";

				ostr << ostrVertices.str();
				ostr << "</IndexedFaceSet>\n";

			X3dElem<t_real>::WriteColor(ostr, m_vecColor);
			ostr << "</Shape>\n";
		}


		void AddVertex(const t_vec& vec) { m_vertices.push_back(vec); }
		void SetColor(const t_vec& vecCol) { m_vecColor = vecCol; }
};

// -----------------------------------------------------------------------------

// lines
template<class t_real = double>
class X3dLines : public X3dElem<t_real>
{
	public:
		using typename X3dElem<t_real>::t_vec;
		using typename X3dElem<t_real>::t_mat;
		using typename X3dElem<t_real>::t_quat;

	protected:
		bool m_bCloseLines = 1;
		std::vector<t_vec> m_vertices;
		t_vec m_vecColor;

	public:
		X3dLines() = default;
		virtual ~X3dLines() = default;

		virtual void Write(std::ostream& ostr) const override
		{
			std::ostringstream ostrVertices;
			ostrVertices << "<Coordinate point=\"";
			for(std::size_t iVert=0; iVert<m_vertices.size(); ++iVert)
			{
				const t_vec& vec = m_vertices[iVert];

				for(t_real d : vec)
					ostrVertices << d << " ";

				if(iVert+1 < m_vertices.size())
					ostrVertices << ", ";
			}
			ostrVertices << "\" />\n";


			ostr << "<Shape>\n";
				ostr << "<IndexedLineSet ";
				ostr << "coordIndex=\"";
				for(std::size_t iVert=0; iVert<m_vertices.size(); ++iVert)
					ostr << iVert << ", ";
				if(m_bCloseLines)
					ostr << 0 << ", ";
				ostr << "-1\"";
				ostr << ">\n";

				ostr << ostrVertices.str();
				ostr << "</IndexedLineSet>\n";

			X3dElem<t_real>::WriteColor(ostr, m_vecColor);
			ostr << "</Shape>\n";
		}

		void AddVertex(const t_vec& vec) { m_vertices.push_back(vec); }
		void SetColor(const t_vec& vecCol) { m_vecColor = vecCol; }

		void SetCloseLines(bool b) { m_bCloseLines = b; }
};

// -----------------------------------------------------------------------------


template<class t_real = double>
class X3d
{
	protected:
		X3dScene<t_real> m_scene;
		std::string m_strComment;

	public:
		X3d() = default;
		virtual ~X3d() = default;

		void SetComment(const std::string& str) { m_strComment = str; }

		void Write(std::ostream& ostr) const
		{
			if(m_strComment.size())
				ostr << "<!--\n" << m_strComment << "\n-->\n";

			ostr << "<X3D>\n";
			m_scene.Write(ostr);
			ostr << "</X3D>\n";
		}

		bool Save(const char* pcFile) const
		{
			std::ofstream ofstr(pcFile);
			if(!ofstr.is_open())
				return false;

			Write(ofstr);
			return true;
		}

		X3dScene<t_real>& GetScene() { return m_scene; }
		const X3dScene<t_real>& GetScene() const { return m_scene; }
};

}
#endif
