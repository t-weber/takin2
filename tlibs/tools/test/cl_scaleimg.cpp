/**
 * tlibs test file
 * @author Tobias Weber <tobias.weber@tum.de>
 * @license GPLv2 or GPLv3
 */

// gcc -o cl_scaleimg cl_scaleimg.cpp ../log/log.cpp -std=c++11 -lOpenCL -lstdc++ -ljpeg -lpng -ltiff

#include "../cl/cl.h"
#include "../gfx/gil.h"
#include "../log/debug.h"
#include "../string/string.h"

#include <iostream>
#include <unordered_map>
#include <type_traits>

namespace gil = tl::gil;

int main(int argc, char** argv)
{
	using t_real = float;


	if(argc < 4)
	{
		std::cerr << "Usage:\n" 
			<< argv[0] << " <input image> <output image> <scale factor>"
			<< std::endl;
		return -1;
	}

	std::string strInFile = argv[1];
	std::string strOutFile = argv[2];
	t_real dScale = tl::str_to_var<t_real, std::string>(argv[3]);

	cl::Platform plat;
	cl::Device dev;
	cl::Context ctx;

	if(!tl::get_best_cl_dev<t_real>(plat, dev, &ctx))
	{
		std::cerr << "Cannot get devices." << std::endl;
		return -1;
	}

	tl::log_info("Using: ", tl::trimmed(plat.getInfo<CL_PLATFORM_NAME>()), ", ",
		tl::trimmed(dev.getInfo<CL_DEVICE_NAME>()));


	std::string strSrc = tl::get_cl_typedefs<t_real>() +
	R"RAWSTR(

	__constant sampler_t sam = CLK_NORMALIZED_COORDS_TRUE | CLK_FILTER_NEAREST | CLK_ADDRESS_CLAMP;

	__kernel void tst_img(__read_only image2d_t imgIn, __write_only image2d_t imgOut)
	{
		const int2 iXY = { get_global_id(0), get_global_id(1) };
		const int2 iDimOut = get_image_dim(imgOut);

		const t_real2 fIn = { (t_real)(iXY[0])/((t_real)(iDimOut[0])),
			(t_real)(iXY[1])/((t_real)(iDimOut[1])) };

		//const t_real4 fPix = read_imagef(imgIn, sam, fIn);
		//write_imagef(imgOut, iXY, fPix);
		const uint4 iPix = read_imageui(imgIn, sam, fIn);
		write_imageui(imgOut, iXY, iPix);
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


	auto pimg = tl::load_image<gil::rgb8_image_t>(strInFile.c_str());
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

	const unsigned int iWOut = iW*dScale;
	const unsigned int iHOut = iH*dScale;

	tl::log_info(strInFile, " (", iW, " x ", iH, ") -> ",
		strOutFile, " (", iWOut, " x ", iHOut, ")");


	bool bInOk, bOutOk;
	cl::Image2D imgIn, imgOut;
	bInOk = tl::create_cl_readimg(ctx, vecImg, imgIn, iW, iH, cl::ImageFormat(CL_RGBA, CL_UNSIGNED_INT8));
	bOutOk = tl::create_cl_writeimg(ctx, imgOut, iWOut, iHOut, cl::ImageFormat(CL_RGBA, CL_UNSIGNED_INT8));
	if(!bInOk || !bOutOk)
	{
		std::cerr << "Could not create input or output buffer." << std::endl;
		return -1;
	}

	const std::size_t iPitch = imgOut.getImageInfo<CL_IMAGE_ROW_PITCH>();
	std::vector<unsigned char> vecOut;
	vecOut.resize(iWOut*iHOut*iC);


	kern.setArg(0, imgIn);
	kern.setArg(1, imgOut);

	cl::CommandQueue qu(ctx, dev);

	cl::Event evtKern;
	cl::NDRange ranOffs, ranGlob(iWOut, iHOut), ranLoc;
	qu.enqueueNDRangeKernel(kern, ranOffs, ranGlob, ranLoc, 0, &evtKern);

	cl::Event evtOut;
	std::vector<cl::Event> vecEvts = {evtKern};
	cl::size_t<3> arrOrigin, arrRange;
	arrRange[0] = iWOut; arrRange[1] = iHOut; arrRange[2] = 1;
	qu.enqueueReadImage(imgOut, false, arrOrigin, arrRange, 0, 0, vecOut.data(), &vecEvts, &evtOut);
	evtOut.wait();

	gil::rgba8_view_t _viewOut = gil::interleaved_view(iWOut, iHOut, (gil::rgba8_pixel_t*)vecOut.data(), iPitch);
	auto viewOut = gil::color_converted_view<gil::rgb8_pixel_t>(_viewOut);
	tl::save_view(strOutFile.c_str(), &viewOut);

	return 0;
}
