/**
 * @author Tobias Weber <tobias.weber@tum.de>
 * @license GPLv2
 *
 * ----------------------------------------------------------------------------
 * Takin (inelastic neutron scattering software package)
 * Copyright (C) 2017-2021  Tobias WEBER (Institut Laue-Langevin (ILL),
 *                          Grenoble, France).
 * Copyright (C) 2013-2017  Tobias WEBER (Technische Universitaet Muenchen
 *                          (TUM), Garching, Germany).
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; version 2 of the License.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along
 * with this program; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 * ----------------------------------------------------------------------------
 */

#include <iostream>
#include <vector>
#include <map>
#include "tlibs/math/linalg.h"
#include "tlibs/helper/traits.h"

using namespace tl;

int main()
{
	std::cout << "vec: " << get_type_dim<std::vector<double>>::value << std::endl;
	std::cout << "ublas vec: " << get_type_dim<ublas::vector<double>>::value << std::endl;
	std::cout << "ublas mat: " << get_type_dim<ublas::matrix<double>>::value << std::endl;
	std::cout << "map: " << get_type_dim<std::map<int, double>>::value << std::endl;


	std::cout << "vec: " << int(get_linalg_type<ublas::vector<double>>::value) << std::endl;
	std::cout << "mat: " << int(get_linalg_type<ublas::matrix<double>>::value) << std::endl;
	return 0;
}
