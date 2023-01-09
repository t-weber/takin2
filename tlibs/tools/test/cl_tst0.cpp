/**
 * tlibs test file
 * @author Tobias Weber <tobias.weber@tum.de>
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

// gcc -o cl_tst0 cl_tst0.cpp -std=c++11 -lOpenCL -lstdc++

#include "../cl/cl.h"
#include <iostream>
#include <unordered_map>
#include <type_traits>

int main()
{
	using t_real = double;

	cl::Platform plat;
	cl::Device dev;
	cl::Context ctx;

	if(!tl::get_best_cl_dev<t_real>(plat, dev, &ctx))
	{
		std::cerr << "Cannot get devices." << std::endl;
		return -1;
	}

	std::string strSrc = tl::get_cl_typedefs<t_real>() +
	R"RAWSTR(

	__kernel void tst_angle(__global const t_real4* vec0, __global const t_real4* vec1,
		__global t_real* dRes)
	{
		const t_real len0 = sqrt(dot(*vec0, *vec0));
		const t_real len1 = sqrt(dot(*vec1, *vec1));

		*dRes = acos(dot(*vec0, *vec1) / (len0*len1)) / M_PI*180.;
	}

	)RAWSTR";


	cl::Program prog;
	std::unordered_map<std::string, cl::Kernel> mapKerns;
	if(!tl::build_cl_src<t_real>(ctx, strSrc, prog, mapKerns))
	{
		std::cerr << "Cannot build kernel." << std::endl;
		return -1;
	}


	cl::Kernel& kern = mapKerns["tst_angle"];

	std::vector<t_real> vecIn = {1., 2., 3., 4.};
	cl::Buffer bufIn(ctx, CL_MEM_READ_ONLY|CL_MEM_COPY_HOST_PTR,
		vecIn.size()*sizeof(decltype(vecIn)::value_type), vecIn.data());
	std::vector<t_real> vecIn2 = {2., 2., 4., 4.};
	cl::Buffer bufIn2(ctx, CL_MEM_READ_ONLY|CL_MEM_COPY_HOST_PTR,
		vecIn2.size()*sizeof(decltype(vecIn2)::value_type), vecIn2.data());
	cl::Buffer bufOut(ctx, CL_MEM_WRITE_ONLY, sizeof(t_real), nullptr);


	kern.setArg(0, bufIn);
	kern.setArg(1, bufIn2);
	kern.setArg(2, bufOut);


	cl::CommandQueue qu(ctx, dev);

	cl::Event evtKern, evtOut;
	cl::NDRange ranOffs, ranGlob(4), ranLoc;
	qu.enqueueNDRangeKernel(kern, ranOffs, ranGlob, ranLoc, 0, &evtKern);

	t_real dOut;
	std::vector<cl::Event> vecEvts = {evtKern};
	qu.enqueueReadBuffer(bufOut, 0, 0, sizeof dOut, &dOut, &vecEvts, &evtOut);
	evtOut.wait();

	std::cout << dOut << std::endl;
	return 0;
}
