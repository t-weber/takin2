/**
 * tlibs test file
 * @author Tobias Weber <tobias.weber@tum.de>
 * @license GPLv2 or GPLv3
 */

// gcc -o cl_tst2 cl_tst2.cpp ../log/log.cpp -std=c++11 -lOpenCL -lstdc++ -ljpeg -lpng -ltiff

#include "../cl/cl.h"
#include "../gfx/gil.h"
#include "../log/debug.h"

#include <iostream>
#include <unordered_map>
#include <type_traits>

namespace gil = tl::gil;

int main()
{
	using t_real = float;

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

	__kernel void tst_img_rgba(uint iW, uint iH,
		__global const uchar4 *pimgIn, __global uchar4 *pimgOut)
	{
		const uint2 iXY = { get_global_id(0), get_global_id(1) };

		uchar4 iPix = pimgIn[iXY.y*iW + iXY.x];
		pimgOut[iXY.y*iW + iXY.x] = iPix.xxxw;
	}

	)RAWSTR";


	cl::Program prog;
	std::unordered_map<std::string, cl::Kernel> mapKerns;
	if(!tl::build_cl_src<t_real>(ctx, strSrc, prog, mapKerns))
	{
		std::cerr << "Cannot build kernel." << std::endl;
		return -1;
	}


	cl::Kernel& kern = mapKerns["tst_img_rgba"];


	auto pimg = tl::load_image<gil::rgb8_image_t>("/home/tweber/Pictures/rac.jpg");
	if(!pimg)
	{
		std::cerr << "Cannot load image." << std::endl;
		return -1;
	}
	auto view = gil::color_converted_view<gil::rgba8_pixel_t>(gil::view(*pimg));
	const unsigned int iW = view.width();
	const unsigned int iH = view.height();
	const unsigned int iC = view.num_channels();
	auto vecImg = tl::view_to_container(view);

	std::cout << "image: " << iW << " x " << iH << " x " << iC << std::endl;

	std::vector<unsigned char> vecOut;
	vecOut.resize(iW*iH*iC);

	cl::Buffer imgIn, imgOut;
	bool bInOk = tl::create_cl_readbuf(ctx, vecImg, imgIn);
	bool bOutOk = tl::create_cl_writebuf<decltype(vecOut)::value_type>(ctx, vecOut.size(), imgOut);

	if(!bInOk || !bOutOk)
	{
		std::cerr << "Could not create input or output buffer." << std::endl;
		return -1;
	}


	kern.setArg(0, iW);
	kern.setArg(1, iH);
	kern.setArg(2, imgIn);
	kern.setArg(3, imgOut);

	cl::CommandQueue qu(ctx, dev);

	cl::Event evtKern;
	cl::NDRange ranOffs, ranGlob(iW, iH), ranLoc;
	std::cout << "Running kernel." << std::endl;
	qu.enqueueNDRangeKernel(kern, ranOffs, ranGlob, ranLoc, 0, &evtKern);

	cl::Event evtOut;
	std::vector<cl::Event> vecEvts = {evtKern};
	cl::size_t<3> arrOrigin, arrRange;
	arrRange[0] = iW; arrRange[1] = iH; arrRange[2] = 1;
	qu.enqueueReadBuffer(imgOut, 0, 0, vecOut.size()*sizeof(decltype(vecOut)::value_type), vecOut.data(), &vecEvts, &evtOut);
	std::cout << "Reading buffer." << std::endl;
	evtOut.wait();

	gil::rgba8_view_t _viewOut = gil::interleaved_view(iW, iH, (gil::rgba8_pixel_t*)vecOut.data(), iW*iC);
	auto viewOut = gil::color_converted_view<gil::rgb8_pixel_t>(_viewOut);
	tl::save_view("/tmp/tst2.png", &viewOut);

	return 0;
}
