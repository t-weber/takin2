/**
 * test dynamic loading of Julia
 *
 * @author Tobias Weber <tweber@ill.fr>
 * @date 17/feb/2020
 * @license GPLv2 or GPLv3
 *
 * g++ -I/usr/local/include/julia -o tst_dl tst_dl.cpp -lboost_filesystem -lboost_system -ldl
 */


#include "../jl.h"
#include <iostream>


int main()
{
	tl::LibJulia jl;

	jl.jl_init();
	std::cout << jl.jl_ver_string() << std::endl;

	jl.jl_eval_string("println(123*100);");
	return 0;
}
