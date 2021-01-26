/**
 * simplified x3d file handling
 * @author Tobias Weber <tobias.weber@tum.de>
 * @date nov-2015
 * @license GPLv2 or GPLv3
 */

#include "x3d.h"


namespace tl{
	template class X3dElem<double>;
	template class X3dScene<double>;
	template class X3dTrafo<double>;
	template class X3dSphere<double>;
	template class X3dCube<double>;
	template class X3dCylinder<double>;
	template class X3dPolygon<double>;
	template class X3dLines<double>;
	template class X3d<double>;
}
