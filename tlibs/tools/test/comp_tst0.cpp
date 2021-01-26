/**
 * boost.comp test
 * @author Tobias Weber <tobias.weber@tum.de>
 * @date jun-2017
 * @license GPLv2 or GPLv3
 */

// gcc -o comp_tst0 comp_tst0.cpp -std=c++11 -lOpenCL -lstdc++ -lm

#include "../cl/comp.h"
#include <iostream>
#include <boost/compute/algorithm/accumulate.hpp>
#include <boost/compute/algorithm/transform.hpp>
#include <boost/compute/algorithm/remove_if.hpp>
#include "../math/math.h"


using namespace tl;
using t_real = double;
using t_vec4 = t_comp_vec4<t_real>;


int main()
{
	// ------------------------------------------------------------------------
	// iterate devices
	std::vector<comp::device> vecDevs = comp::system::devices();
	std::cout << "Found " << vecDevs.size() << " CL devices:" << std::endl;

	for(const comp::device& dev : vecDevs)
	{
		std::cout << "\t * ";
		if(dev.type() & comp::device::gpu) std::cout << "GPU";
		if(dev.type() & comp::device::cpu) std::cout << "CPU";
		std::cout << ": " << dev.global_memory_size()/1024/1024 << " MiB";
		std::cout << ", " << dev.clock_frequency() << " MHz";
		std::cout << ", " << dev.compute_units() << " units";
		std::cout << " (" << dev.name() << ")" << std::endl;
	}
	std::cout << std::endl;
	// ------------------------------------------------------------------------



	// ------------------------------------------------------------------------
	// get device
	comp::device dev = *tl::get_best_comp_dev<t_real>();
	std::cout << "Using " << dev.name() << ".\n" << std::endl;

	// context & queue
	comp::context cont(dev);
	comp::command_queue qu(cont, dev);
	// ------------------------------------------------------------------------



	// ------------------------------------------------------------------------
	// data
	const std::size_t iSize = 64;
	std::vector<t_real> vec = tl::linspace<t_real>(0., 100., iSize);

	std::cout << "in: ";
	for(t_real d : vec) std::cout << d << ", ";
	std::cout << std::endl;
	// ------------------------------------------------------------------------



	// ------------------------------------------------------------------------
	// copy to GPU
	comp::vector<t_real> vecCL(iSize, cont);
	auto futCpy1 = comp::copy_async(vec.begin(), vec.end(), vecCL.begin(), qu);
	std::cout << "Copying data to GPU ... ";
	futCpy1.wait();
	std::cout << "OK." << std::endl << std::endl;
	// ------------------------------------------------------------------------




	// ------------------------------------------------------------------------
	// testing higher-order function 1: accumulate
	BOOST_COMPUTE_FUNCTION(t_real, fkt_accum, (t_real d1, t_real d2), { return d1 + 2.*d2; } );

	t_real dSum1 = std::accumulate(vec.begin(), vec.end(), t_real(0.),
		[](t_real d1, t_real d2) -> t_real { return d1+2.*d2; });
	t_real dSum2 = comp::accumulate(vecCL.begin(), vecCL.end(), t_real(0.),
		fkt_accum, qu);
	std::cout << "accum (CPU): " << dSum1 << std::endl;
	std::cout << "accum (GPU): " << dSum2 << std::endl;
	std::cout << std::endl;
	// ------------------------------------------------------------------------



	// ------------------------------------------------------------------------
	// testing higher-order function 2: transform
	BOOST_COMPUTE_FUNCTION(t_real, fkt_trafo, (t_real d), { return sqrt(d) + 10.; } );

	std::vector<t_real> vecTst(iSize);
	std::transform(vec.begin(), vec.end(), vecTst.begin(),
		[](t_real d) -> t_real { return std::sqrt(d)+10.; } );
	comp::transform(vecCL.begin(), vecCL.end(), vecCL.begin(), fkt_trafo, qu);

	auto futCpy2 = comp::copy_async(vecCL.begin(), vecCL.end(), vec.begin(), qu);
	std::cout << "Fetching data from GPU ... ";
	futCpy2.wait();
	std::cout << "OK." << std::endl;

	std::cout << "trafo (CPU): ";
	for(t_real d : vecTst) std::cout << d << ", ";
	std::cout << std::endl;

	std::cout << "trafo (GPU): ";
	for(t_real d : vec) std::cout << d << ", ";
	std::cout << std::endl << std::endl;
	// ------------------------------------------------------------------------



	// ------------------------------------------------------------------------
	// testing higher-order function 3: remove_if
	BOOST_COMPUTE_FUNCTION(bool, fkt_rem, (t_real d), { return d > 15.; } );

	auto iterEnd = std::remove_if(vec.begin(), vec.end(),
		[](t_real d) -> bool { return d > 15.; } );
	auto iterEndGPU = comp::remove_if(vecCL.begin(), vecCL.end(), fkt_rem, qu);

	std::cout << "remove_if (CPU): ";
	for(auto iter = vec.begin(); iter!=iterEnd; ++iter) std::cout << *iter << ", ";
	std::cout << std::endl;

	auto futCpy3 = comp::copy_async(vecCL.begin(), iterEndGPU, vec.begin(), qu);
	std::cout << "Fetching data from GPU ... ";
	futCpy3.wait();
	std::cout << "OK." << std::endl;
	vec.resize(iterEndGPU - vecCL.begin());

	std::cout << "remove_if (GPU): ";
	for(t_real d : vec) std::cout << d << ", ";
	std::cout << std::endl << std::endl;
	// ------------------------------------------------------------------------



	// ------------------------------------------------------------------------
	// 4-vector test
	BOOST_COMPUTE_FUNCTION(t_vec4, fkt_accum_vec, (const t_vec4& vec1, const t_vec4& vec2),
	{
		const t_real4 vec3 = {0,0,1,0};
		return vec1 + cross(vec2, vec3);
	});

	inject_comp_typedefs<t_real>(fkt_accum_vec);
	//std::cout << fkt_accum_vec.source() << std::endl;


	std::vector<std::vector<t_real>> vecs4 = 
	{
		{{1., 2., 3., 0.}},
		{{1., 2., 3., 0.}},
		{{1., 2., 3., 0.}},
		{{1., 2., 3., 0.}},
	};

	comp::vector<t_vec4> vecs4CL(vecs4.size(), cont);
	std::vector<comp::future<void>> vecFuts;

	for(std::size_t iIdx=0; iIdx<vecs4.size(); ++iIdx)
	{
		auto fut = comp::copy_async((t_vec4*)vecs4[iIdx].data(),
			(t_vec4*)(vecs4[iIdx].data() + vecs4[iIdx].size()),
			vecs4CL.begin()+iIdx, qu);
		vecFuts.emplace_back(std::move(fut));
	}
	for(std::size_t iIdx=0; iIdx<vecFuts.size(); ++iIdx)
		vecFuts[iIdx].wait();

	t_vec4 vecSum = comp::accumulate(vecs4CL.begin(), vecs4CL.end(),
		t_vec4(0,0,0,0), fkt_accum_vec, qu);
	std::cout << "vec accum: " << vecSum << std::endl;
	// ------------------------------------------------------------------------


	return 0;
}
