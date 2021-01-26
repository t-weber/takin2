/**
 * magnetic space group library
 * @author Tobias Weber <tweber@ill.fr>
 * @date 7-apr-18
 * @license GPLv3, see 'LICENSE' file
 * @desc The present version was forked on 8-Nov-2018 from my privately developed "magtools" project (https://github.com/t-weber/magtools).
 */

#include "libs/magsg.h"
#include "libs/math20.h"
using namespace m_ops;


int main()
{
	using t_real = double;
	using t_vec = std::vector<t_real>;
	using t_mat = tl2::mat<t_real, std::vector>;

	Spacegroups<t_mat, t_vec> sgs;
	sgs.Load("magsg.xml");

	return 0;
}
