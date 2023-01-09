/**
 * tlibs2 -- swig interface
 * @author Tobias Weber <tweber@ill.fr>
 * @date 4-jun-2020
 * @license see 'LICENSE' file
 *
 * ----------------------------------------------------------------------------
 * tlibs
 * Copyright (C) 2017-2021  Tobias WEBER (Institut Laue-Langevin (ILL),
 *                          Grenoble, France).
 * Copyright (C) 2015-2017  Tobias WEBER (Technische Universitaet Muenchen
 *                          (TUM), Garching, Germany).
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, version 3 of the License.
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

%module tl2
%{
	//#include "magdyn.h"
	#include "instr.h"
%}

%include "std_vector.i"
%include "std_array.i"
%include "std_map.i"
%include "std_unordered_map.i"
%include "std_pair.i"
//%include "std_tuple.i"
%include "std_string.i"

%template(VecStr) std::vector<std::string>;
%template(VecD) std::vector<double>;
%template(ArrD3) std::array<double,3>;
%template(ArrD5) std::array<double,5>;
%template(ArrB3) std::array<bool,3>;
%template(MapStrStr) std::unordered_map<std::string, std::string>;

//%include "magdyn.h"
%include "instr.h"

%template(FileInstrBaseD) tl2::FileInstrBase<double>;
