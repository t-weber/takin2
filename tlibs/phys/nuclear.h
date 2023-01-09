/**
 * nuclear physics formulas
 * @author Tobias Weber <tobias.weber@tum.de>
 * @date 2016
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

#ifndef __TLIBS_NUCLEAR__
#define __TLIBS_NUCLEAR__


#include "units.h"
#include "../math/math.h"
#include "../helper/exception.h"
#include <cmath>


namespace tl {

/**
 * see: https://de.wikipedia.org/wiki/Mittleres_logarithmisches_Energiedekrement
 */
template<class Y=double>
Y mean_log_E_loss(Y A)
{
	Y a = Y((1-A)*(1-A)) / Y((1+A)*(1+A));
	return Y(1) + a/(Y(1)-a) * std::log(a);
}


/**
 * see: https://de.wikipedia.org/wiki/Mittleres_logarithmisches_Energiedekrement
 */
template<class Sys, class Y=double>
Y mean_collisions(Y A, const t_energy<Sys,Y>& E_from, const t_energy<Sys,Y>& E_to)
{
	Y dLogLoss = mean_log_E_loss<Y>(A);
	Y dLogE = std::log(Y(E_from/E_to));

	return dLogE / dLogLoss;
}

}

#endif
