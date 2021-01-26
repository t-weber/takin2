/**
 * Form factor and scattering length tables
 * @author Tobias Weber <tobias.weber@tum.de>
 * @date nov-2015
 * @license GPLv2
 */

#include "formfact.h"
#include "formfact_impl.h"

namespace xtl {

template class PeriodicSystem<double>;
template class FormfactList<double>;
template class MagFormfactList<double>;
template class ScatlenList<double>;

template class PeriodicSystem<float>;
template class FormfactList<float>;
template class MagFormfactList<float>;
template class ScatlenList<float>;

}
