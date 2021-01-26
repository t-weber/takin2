/**
 * Loads tabulated spacegroups
 * @author Tobias Weber <tobias.weber@tum.de>
 * @date feb-2016
 * @license GPLv2
 */

#define NO_QT
#include "spacegroup.h"
#include "spacegroup_impl.h"


namespace xtl {

template class SpaceGroup<double>;
template class SpaceGroups<double>;

template class SpaceGroup<float>;
template class SpaceGroups<float>;

}
