/**
 * tlibs test file
 * @author Tobias Weber <tobias.weber@tum.de>
 * @license GPLv2 or GPLv3
 */

// gcc -o cl_tst1 cl_tst1.cpp ../log/log.cpp -std=c++11 -lOpenCL -lstdc++ -ljpeg -lpng -ltiff

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

	std::cout << "Using: " << plat.getInfo<CL_PLATFORM_NAME>() << ", "
		<< dev.getInfo<CL_DEVICE_NAME>() << std::endl;


	std::string strSrc = tl::get_cl_typedefs<t_real>() +
	R"RAWSTR(

	__constant sampler_t sam = CLK_NORMALIZED_COORDS_FALSE | CLK_FILTER_NEAREST | CLK_ADDRESS_CLAMP;
	//__constant sampler_t sam = CLK_NORMALIZED_COORDS_TRUE | CLK_FILTER_LINEAR | CLK_ADDRESS_CLAMP;

	__kernel void tst_img(__read_only image2d_t imgIn, __write_only image2d_t imgOut)
	{
		const int2 iXY = { get_global_id(0), get_global_id(1) };
		const uint4 iPix = read_imageui(imgIn, sam, iXY);

		uint4 iPixOut = {iPix[0]/2, iPix[1]/2, iPix[2]/2, 0xff};
		write_imageui(imgOut, iXY, iPixOut);

		//write_imagef(imgOut, iXY, (t_real4)(1.,1.,1.,1.));
	}

	)RAWSTR";


	cl::Program prog;
	std::unordered_map<std::string, cl::Kernel> mapKerns;
	if(!tl::build_cl_src<t_real>(ctx, strSrc, prog, mapKerns))
	{
		std::cerr << "Cannot build kernel." << std::endl;
		return -1;
	}


	cl::Kernel& kern = mapKerns["tst_img"];


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

/*	int iW = 2, iH = 2;
	std::vector<unsigned char> vecImg =
		{0xff,0xff,0xff,0xff,
		0x00,0x00,0x00,0xff,
		0x00,0x00,0x00,0xff,
		0xff,0xff,0xff,0xff};*/

	std::cout << "image: " << iW << " x " << iH << " x " << iC << std::endl;
	//std::cout << tl::get_typename<decltype(vecImg)>() << std::endl;
	//std::cout << vecImg.size() << std::endl;
	//std::cout << iW*iH << std::endl;
	//std::cout << sizeof(decltype(vecImg)::value_type) << std::endl;
	//std::cout << "data type is pod: " << std::is_pod<decltype(vecImg)::value_type>::value << std::endl;

	bool bInOk, bOutOk;
	cl::Image2D imgIn, imgOut;
	bInOk = tl::create_cl_readimg(ctx, vecImg, imgIn, iW, iH, cl::ImageFormat(CL_RGBA, CL_UNSIGNED_INT8));
	bOutOk = tl::create_cl_writeimg(ctx, imgOut, iW, iH, cl::ImageFormat(CL_RGBA, CL_UNSIGNED_INT8));
	if(!bInOk || !bOutOk)
	{
		std::cerr << "Could not create input or output buffer." << std::endl;
		return -1;
	}

	//std::cout << imgIn.getImageInfo<CL_IMAGE_WIDTH>() << " " << imgIn.getImageInfo<CL_IMAGE_HEIGHT>() << std::endl;
	const std::size_t iPitch = imgOut.getImageInfo<CL_IMAGE_ROW_PITCH>();
	std::cout << "row pitch: " << iPitch << " bytes" << std::endl;


	std::vector<unsigned char> vecOut;
	vecOut.resize(iW*iH*iC);


	kern.setArg(0, imgIn);
	kern.setArg(1, imgOut);

	cl::CommandQueue qu(ctx, dev);

	cl::Event evtKern;
	cl::NDRange ranOffs, ranGlob(iW, iH), ranLoc;
	std::cout << "Running kernel." << std::endl;
	qu.enqueueNDRangeKernel(kern, ranOffs, ranGlob, ranLoc, 0, &evtKern);

	cl::Event evtOut;
	std::vector<cl::Event> vecEvts = {evtKern};
	cl::size_t<3> arrOrigin, arrRange;
	arrRange[0] = iW; arrRange[1] = iH; arrRange[2] = 1;
	qu.enqueueReadImage(imgOut, false, arrOrigin, arrRange, 0, 0, vecOut.data(), &vecEvts, &evtOut);
	std::cout << "Reading buffer." << std::endl;
	evtOut.wait();

	//std::fill(vecOut.begin(), vecOut.end(), 0xff);
	gil::rgba8_view_t _viewOut = gil::interleaved_view(iW, iH, (gil::rgba8_pixel_t*)vecOut.data(), iPitch);
	auto viewOut = gil::color_converted_view<gil::rgb8_pixel_t>(_viewOut);
	tl::save_view("/tmp/tst1.png", &viewOut);

	return 0;
}
