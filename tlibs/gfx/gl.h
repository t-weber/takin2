/**
 * GL drawing
 * @author Tobias Weber <tobias.weber@tum.de>
 * @date 22-dec-2014
 * @license GPLv2 or GPLv3
 */

#ifndef __TLIBS_GL_STUFF_H__
#define __TLIBS_GL_STUFF_H__

#include <string>
#include <ostream>
#include <unordered_map>
#include <gl.h>
//#include <glext.h>

#include "../math/linalg.h"
#include "../math/geo.h"

namespace tl {

// ----------------------------------------------------------------------------
/**
 * GL type traits
 */
template<class T = GLdouble> struct gl_traits {};

template<> struct gl_traits<GLdouble>
{
	using value_type = GLdouble;

	static void GetProjMatrix(value_type *pdMat) { glGetDoublev(GL_PROJECTION_MATRIX, pdMat); }
	static void GetModelMatrix(value_type *pdMat) { glGetDoublev(GL_MODELVIEW_MATRIX, pdMat); }

	static void LoadMatrix(value_type* pdMat) { glLoadMatrixd(pdMat); }
	static void MultMatrix(value_type* pdMat) { glMultMatrixd(pdMat); }

	static void SetScale(value_type x, value_type y, value_type z) { glScaled(x, y, z); }
	static void SetTranslate(value_type x, value_type y, value_type z) { glTranslated(x, y, z); }

	static void SetVertex(value_type x, value_type y) { glVertex2d(x, y); }
	static void SetVertex(value_type x, value_type y, value_type z) { glVertex3d(x, y, z); }
	static void SetVertex(const value_type *px) { glVertex3dv(px); }
	static void SetNorm(value_type x, value_type y, value_type z) { glNormal3d(x, y, z); }
	static void SetNorm(const value_type* px) { glNormal3dv(px); }
	static void SetTextureCoord(value_type u, value_type v) { glTexCoord2d(u, v); }

	static void SetColor(value_type r, value_type g, value_type b) { glColor3d(r, g, b); }
	static void SetColor(value_type r, value_type g, value_type b, value_type a) { glColor4d(r, g, b, a); }
	static void SetMaterial(GLenum param1, GLenum param2, const value_type* pVal)
	{ // not available as "double"
		using T = GLfloat;
		T tconv[] = { T(pVal[0]), T(pVal[1]), T(pVal[2]), T(pVal[3]) };
		glMaterialfv(param1, param2, tconv);
	}
	static void SetLight(GLenum param1, GLenum param2, const value_type *pVec)
	{ // not available as "double"
		using T = GLfloat;
		T tconv[] = { T(pVec[0]), T(pVec[1]), T(pVec[2]), T(pVec[3]) };
		glLightfv(param1, param2, tconv);
	}
};

template<> struct gl_traits<GLfloat>
{
	using value_type = GLfloat;

	static void GetProjMatrix(value_type *pdMat) { glGetFloatv(GL_PROJECTION_MATRIX, pdMat); }
	static void GetModelMatrix(value_type *pdMat) { glGetFloatv(GL_MODELVIEW_MATRIX, pdMat); }

	static void LoadMatrix(value_type* pdMat) { glLoadMatrixf(pdMat); }
	static void MultMatrix(value_type* pdMat) { glMultMatrixf(pdMat); }

	static void SetScale(value_type x, value_type y, value_type z) { glScalef(x, y, z); }
	static void SetTranslate(value_type x, value_type y, value_type z) { glTranslatef(x, y, z); }

	static void SetVertex(value_type x, value_type y) { glVertex2f(x, y); }
	static void SetVertex(value_type x, value_type y, value_type z) { glVertex3f(x, y, z); }
	static void SetVertex(const value_type *px) { glVertex3fv(px); }
	static void SetNorm(value_type x, value_type y, value_type z) { glNormal3f(x, y, z); }
	static void SetNorm(const value_type* px) { glNormal3fv(px); }
	static void SetTextureCoord(value_type u, value_type v) { glTexCoord2f(u, v); }

	static void SetColor(value_type r, value_type g, value_type b) { glColor3f(r, g, b); }
	static void SetColor(value_type r, value_type g, value_type b, value_type a) { glColor4f(r, g, b, a); }
	static void SetMaterial(GLenum param1, GLenum param2, const value_type* pVal) { glMaterialfv(param1, param2, pVal); }
	static void SetLight(GLenum param1, GLenum param2, const value_type *pVec) { glLightfv(param1, param2, pVec); }
};


template<class T = GLdouble>
using t_mat4_gen = ublas::matrix<T, ublas::row_major, ublas::bounded_array<T,4*4>>;
template<class T = GLdouble>
using t_mat3_gen = ublas::matrix<T, ublas::row_major, ublas::bounded_array<T,3*3>>;
template<class T = GLdouble>
using t_vec4_gen = ublas::vector<T, ublas::bounded_array<T,4>>;
template<class T = GLdouble>
using t_vec3_gen = ublas::vector<T, ublas::bounded_array<T,3>>;

typedef t_mat4_gen<GLdouble> t_mat4d;
typedef t_mat3_gen<GLdouble> t_mat3d;
typedef t_vec4_gen<GLdouble> t_vec4d;
typedef t_vec3_gen<GLdouble> t_vec3d;

typedef t_mat4_gen<GLfloat> t_mat4f;
typedef t_mat3_gen<GLfloat> t_mat3f;
typedef t_vec4_gen<GLfloat> t_vec4f;
typedef t_vec3_gen<GLfloat> t_vec3f;
// ----------------------------------------------------------------------------



template<typename T=GLdouble, typename... Args>
void to_gl_array(const ublas::matrix<T, Args...>& mat, T* glmat)
{
	glmat[0]=mat(0,0);  glmat[1]=mat(1,0);  glmat[2]=mat(2,0);
	glmat[4]=mat(0,1);  glmat[5]=mat(1,1);  glmat[6]=mat(2,1);
	glmat[8]=mat(0,2);  glmat[9]=mat(1,2);  glmat[10]=mat(2,2);

	if(mat.size1()>=4 && mat.size2()>=4)
	{
		glmat[3]=mat(3,0); glmat[7]=mat(3,1); glmat[11]=mat(3,2);
		glmat[12]=mat(0,3); glmat[13]=mat(1,3); glmat[14]=mat(2,3); glmat[15]=mat(3,3);
	}
	else
	{
		glmat[3]=0; glmat[7]=0; glmat[11]=0;
		glmat[12]=0; glmat[13]=0; glmat[14]=0; glmat[15]=1;
	}
}

template<typename t_mat = t_mat4_gen<GLdouble>>
t_mat from_gl_array(const typename t_mat::value_type* glmat)
{
	t_mat mat(4,4);

	for(short j=0; j<4; ++j)
		for(short i=0; i<4; ++i)
			mat(i,j)=glmat[i + j*4];

	return mat;
}


/**
 * project a point using given projection and modelview matrices
 */
template<typename t_mat = t_mat4_gen<GLdouble>, typename t_vec = t_vec4_gen<GLdouble>,
	typename T = typename t_mat::value_type>
void proj_pt(T dX, T dY, T dZ, const t_mat& matProj, const t_mat& matMV,
	T& dXProj, T& dYProj)
{
	t_mat mat = prod_mm(matProj, matMV);
	t_vec vec = prod_mv(mat, make_vec<t_vec>({dX, dY, dZ, 1.}));
	vec /= vec[3];

	dXProj = vec[0];
	dYProj = vec[1];
}


/**
 * project a point using current projection and modelview matrices
 */
template<typename t_mat = t_mat4_gen<GLdouble>, typename t_vec = t_vec4_gen<GLdouble>,
	typename T = typename t_mat::value_type>
void gl_proj_pt(T dX, T dY, T dZ, T& dXProj, T& dYProj)
{
	T dMatMV[16], dMatProj[16];
	gl_traits<T>::GetProjMatrix(dMatProj);
	gl_traits<T>::GetModelMatrix(dMatMV);

	t_mat matProj = from_gl_array<t_mat>(dMatProj);
	t_mat matMV = from_gl_array<t_mat>(dMatMV);

	proj_pt<t_mat, t_vec, T>(dX, dY, dZ, matProj, matMV, dXProj, dYProj);
}


/**
 * multiply a point with the inverse modelview matrix
 */
template<typename t_mat = t_mat4_gen<GLdouble>, typename t_vec = t_vec4_gen<GLdouble>,
	typename T = typename t_mat::value_type>
void gl_mv_pt(const t_vec& vec, t_vec& vecOut)
{
	T dMatMV[16];
	gl_traits<T>::GetModelMatrix(dMatMV);

	t_mat matMV = from_gl_array<t_mat>(dMatMV);
	t_mat matMV_inv;
	tl::inverse(matMV, matMV_inv);

	vecOut = prod_mv(matMV_inv, vec);
}

/**
 * distance to the object defined by the current modelview matrix
 */
template<typename t_mat = t_mat4_gen<GLdouble>, typename t_vec = t_vec4_gen<GLdouble>,
	typename T = typename t_mat::value_type>
T gl_dist_mv()
{
	t_vec vecPos;
	gl_mv_pt(make_vec<t_vec>({0.,0.,0.,1.}), vecPos);
	vecPos /= vecPos[3];
	vecPos[3] = 0.;
	T dDist = veclen(vecPos);
	return dDist;
}


/**
 * size of projected sphere using fov angle
 */
template<typename t_mat = t_mat4_gen<GLdouble>, typename t_vec = t_vec4_gen<GLdouble>,
	typename T = typename t_mat::value_type>
T gl_proj_sphere_size(T dFOV, T dRadius)
{
	return T(2)*dRadius / (T(2)*gl_dist_mv<t_mat, t_vec, T>() * std::tan(T(0.5)*dFOV));
}

/**
 * size of projected sphere using projection matrix
 */
template<typename t_mat = t_mat4_gen<GLdouble>, typename t_vec = t_vec4_gen<GLdouble>,
	typename T = typename t_mat::value_type>
T gl_proj_sphere_size(T dRadius)
{
	T dMatProj[16];
	gl_traits<T>::GetProjMatrix(dMatProj);
	t_mat matProj = from_gl_array<t_mat>(dMatProj);

	T dDist = gl_dist_mv<t_mat, t_vec, T>();
	t_vec vec1 = make_vec<t_vec>({0., dRadius, dDist, 1.});
	t_vec vec2 = make_vec<t_vec>({0., -dRadius, dDist, 1.});

	t_vec vecProj1 = prod_mv(matProj, vec1); vecProj1 /= vecProj1[3];
	t_vec vecProj2 = prod_mv(matProj, vec2); vecProj2 /= vecProj2[3];
	t_vec vecProj = vecProj2 - vecProj1;

	vecProj[3] = vecProj[2] = 0.;
	return veclen(vecProj);
}


/**
 * ray through screen coordinates (dX, dY)
 * similar to: https://www.opengl.org/sdk/docs/man2/xhtml/gluUnProject.xml
 */
template<typename t_mat = t_mat4_gen<GLdouble>, typename t_vec = t_vec4_gen<GLdouble>,
	typename T = typename t_mat::value_type>
Line<T> screen_ray(T dX, T dY, const t_mat& matProj, const t_mat& matMV)
{
	t_mat mat = prod_mm(matProj, matMV);
	t_mat matInv;
	inverse(mat, matInv);

	t_vec vecNear = make_vec<t_vec>({dX, dY, T(-1.), T(1.)});
	t_vec vecFar = make_vec<t_vec>({dX, dY, T(1.), T(1.)});

	vecNear = prod_mv(matInv, vecNear);
	vecFar = prod_mv(matInv, vecFar);

	vecNear /= vecNear[3];
	vecFar /= vecFar[3];

	ublas::vector<T> vecPos = vecNear; vecPos.resize(3, 1);
	ublas::vector<T> vecDir = vecFar-vecNear; vecDir.resize(3, 1);
	Line<T> line(vecPos, vecDir);
	return line;
}


/**
 * ray through screen coordinates (dX, dY) using current proj & MV matrices
 * similar to: https://www.opengl.org/sdk/docs/man2/xhtml/gluUnProject.xml
 */
template<typename t_mat = t_mat4_gen<GLdouble>, typename t_vec = t_vec4_gen<GLdouble>,
	typename T = typename t_mat::value_type>
Line<T> gl_screen_ray(T dX, T dY)
{
	T dMatMV[16], dMatProj[16];
	gl_traits<T>::GetProjMatrix(dMatProj);
	gl_traits<T>::GetModelMatrix(dMatMV);

	t_mat matProj = from_gl_array<t_mat>(dMatProj);
	t_mat matMV = from_gl_array<t_mat>(dMatMV);

	return screen_ray<t_mat, t_vec, T>(dX, dY, matProj, matMV);
}

// --------------------------------------------------------------------------------


/**
 * simple camera class
 */
template<typename t_mat = ublas::matrix<GLdouble>, typename t_vec = ublas::vector<GLdouble>,
	typename T = typename t_mat::value_type>
class Cam
{
protected:
	t_vec m_vecDir = make_vec<t_vec>({0.,0.,1.});
	t_vec m_vecUp = make_vec<t_vec>({0.,1.,0.});
	t_vec m_vecRight = make_vec<t_vec>({1.,0.,0.});
	t_vec m_vecPos = make_vec<t_vec>({0.,0.,0.});

public:
	Cam() {}
	Cam(const t_vec& vecDir, const t_vec& vecUp)
		: m_vecDir(vecDir), m_vecUp(vecUp)
	{
		m_vecDir = veclen(m_vecDir);
		m_vecUp = veclen(m_vecUp);
		m_vecRight = cross_3(m_vecUp, m_vecDir);
		m_vecUp = cross_3(m_vecDir, m_vecRight);
	}

	virtual ~Cam() {}

	// only rotation matrix
	t_mat GetRotMatrix() const
	{
		return row_matrix({m_vecRight, m_vecUp, m_vecDir});
	}

	// homogeneous coordinates
	template<typename t_mat_h>
	t_mat_h GetHomMatrix() const
	{
		t_mat_h mat = GetRotMatrix();
		mat.resize(4,4,1);
		for(int i=0; i<3; ++i)
			mat(3,i) = mat(i,3) = 0.;
		mat(3,3) = 1.;

		t_mat_h matTrans = make_mat<t_mat_h>(
			{{ 1, 0, 0, m_vecPos[0]},
			{  0, 1, 0, m_vecPos[1]},
			{  0, 0, 1, m_vecPos[2]},
			{  0, 0, 0,          1 }});

		mat = prod_mm(mat, matTrans);
		return mat;
	}

	void RotateAroundRight(T tAngle)
	{
		t_mat matRot = rotation_matrix<t_mat>(m_vecRight, tAngle);
		m_vecUp = prod_vm(m_vecUp, matRot);
		m_vecDir = prod_vm(m_vecDir, matRot);
	}
	void RotateAroundUp(T tAngle)
	{
		t_mat matRot = rotation_matrix<t_mat>(m_vecUp, tAngle);
		m_vecRight = prod_vm(m_vecRight, matRot);
		m_vecDir = prod_vm(m_vecDir, matRot);
	}
	void RotateAroundDir(T tAngle)
	{
		t_mat matRot = rotation_matrix<t_mat>(m_vecDir, tAngle);
		m_vecRight = prod_vm(m_vecRight, matRot);
		m_vecUp = prod_vm(m_vecUp, matRot);
	}

	void MoveForward(T t)
	{
		m_vecPos += m_vecDir * t;
	}

	const t_vec& GetPos() const { return m_vecPos; }
};

// --------------------------------------------------------------------------------
}

//
//#ifdef TLIBS_INC_HDR_IMPLS
//	#include "gl_impl.h"
//#endif

#endif
