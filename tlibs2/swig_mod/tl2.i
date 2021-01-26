/**
 * tlibs2 -- swig interface
 * @author Tobias Weber <tweber@ill.fr>
 * @date 4-jun-2020
 * @license see 'LICENSE' file
 */

%module tl2
%{
	#include "instr.h"
%}

%include "std_vector.i"
%include "std_array.i"
%include "std_map.i"
%include "std_unordered_map.i"
%include "std_pair.i"
%include "std_string.i"

%template(VecStr) std::vector<std::string>;
%template(VecD) std::vector<double>;
%template(ArrD3) std::array<double,3>;
%template(ArrD5) std::array<double,5>;
%template(ArrB3) std::array<bool,3>;
%template(MapStrStr) std::unordered_map<std::string, std::string>;

%include "instr.h"
%template(FileInstrBaseD) tl2::FileInstrBase<double>;
