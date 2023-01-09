/**
 * minimalistic expression evaluator
 * @author Tobias Weber <tobias.weber@tum.de>
 * @date apr-2016
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

#include <boost/integer_fwd.hpp>
#include "eval.h"
#include "eval_impl.h"

namespace tl
{
	template std::pair<bool, double> eval_expr<std::string, double>(const std::string& str) noexcept;
	template std::pair<bool, float> eval_expr<std::string, float>(const std::string& str) noexcept;
	template std::pair<bool, int> eval_expr<std::string, int>(const std::string& str) noexcept;
}
