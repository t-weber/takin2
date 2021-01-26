/**
 * tlibs test file
 * @author Tobias Weber <tobias.weber@tum.de>
 * @license GPLv2 or GPLv3
 * clang -I/usr/include/lapacke -o spin ../test/spin.cpp ../log/log.cpp -lstdc++ -std=c++11 -lm -llapacke -llapack
 */

#define TLIBS_USE_GLOBAL_OPS
#include "../math/linalg.h"
#include "../math/linalg_ops.h"
#include "../math/linalg2.h"
#include "../math/linalg2_impl.h"
#include "../phys/spin.h"

#include <iostream>

using t_real = double;
using namespace tl;
using t_mat = ublas::matrix<std::complex<t_real>>;
using t_vec = ublas::vector<std::complex<t_real>>;

int main()
{
	auto vec = get_spin_matrices();
	auto I = unit_m<t_mat>(2);

	// operator components for spin 1
	t_mat S1x = tensor_prod(vec[0], I);
	t_mat S1y = tensor_prod(vec[1], I);
	t_mat S1z = tensor_prod(vec[2], I);

	// operator components for spin 2
	t_mat S2x = tensor_prod(I, vec[0]);
	t_mat S2y = tensor_prod(I, vec[1]);
	t_mat S2z = tensor_prod(I, vec[2]);

	// two-spin operator
	t_mat S1S2 = S1x*S2x + S1y*S2y + S1z*S2z;

	std::cout << "S1*S2 (direct) = " << S1S2 << std::endl;


	std::vector<t_vec> evecs;
	std::vector<std::complex<t_real>> evals;
	if(!eigenvec_cplx(S1S2, evecs, evals))
		std::cerr << "Cannot calculate eigenvectors." << std::endl;

	for(std::size_t iEV=0; iEV<evecs.size(); ++iEV)
	{
		std::cout << "evec: " << evecs[iEV]
			<< ", eval: " << evals[iEV]
			<< std::endl;
	}
	std::cout << std::endl;


	// ---------------------------------------------------------


	auto vecLadder = get_ladder_ops();
	std::cout << "ladder ops: " << vecLadder << std::endl;
	std::cout << "rot_x 180: " << tl::rot_spin(0, M_PI) << std::endl;
	std::cout << "rot_y 180: " << tl::rot_spin(1, M_PI) << std::endl;
	std::cout << "rot_z 180: " << tl::rot_spin(2, M_PI) << std::endl;
	std::cout << std::endl;

	t_mat C = 0.5 * (tensor_prod(vecLadder[0], I) * tensor_prod(I, vecLadder[1]) +
		tensor_prod(vecLadder[1], I) * tensor_prod(I, vecLadder[0]))
		+ S1z*S2z;
	std::cout << "S1*S2 (via ladder) = " << C << std::endl;
	//std::cout << int(tl::get_linalg_type<decltype(C)>::value) << std::endl;
	//std::cout << int(tl::get_linalg_type<decltype(S1S2)>::value) << std::endl;
	std::cout << std::boolalpha << (C == S1S2) << std::endl;
	std::cout << std::endl;

	// ---------------------------------------------------------


	ublas::matrix<std::complex<t_real>> mRot = tl::rot_spin(2, M_PI/2.);
	ublas::vector<std::complex<t_real>> vecUp(2);
	vecUp[0] = 1; vecUp[1] = 0;

	std::cout << "0 deg: " << vecUp << std::endl;
	std::cout << "90 deg: " << mRot*vecUp << std::endl;
	std::cout << "180 deg: " << mRot*mRot*vecUp << std::endl;
	std::cout << "270 deg: " << mRot*mRot*mRot*vecUp << std::endl;
	std::cout << "360 deg: " << mRot*mRot*mRot*mRot*vecUp << std::endl;
	std::cout << std::endl;

	return 0;
}
