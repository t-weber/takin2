/**
 * Helper functions using Boost.Compute
 * @author Tobias Weber <tobias.weber@tum.de>
 * @date jun-2017
 * @license GPLv2 or GPLv3
 */

#ifndef __CL_COMP_WRAP_H__
#define __CL_COMP_WRAP_H__

#include "cl.h"
#include <vector>

#include <boost/optional.hpp>
#include <boost/compute/system.hpp>
#include <boost/compute/container/vector.hpp>
#include <boost/compute/type_traits/make_vector_type.hpp>

namespace tl {

namespace comp = boost::compute;


template<class t_scalar> 
using t_comp_vec2 = typename comp::make_vector_type<t_scalar, 2>::type;

template<class t_scalar> 
using t_comp_vec4 = typename comp::make_vector_type<t_scalar, 4>::type;


template<typename t_real=double>
boost::optional<comp::device> get_best_comp_dev()
{
	struct _Dev
	{
		comp::device dev;
		cl_device_type devtype;
	};

	std::vector<_Dev> vecAllDevs;

	// enumerate all platforms & devices
	for(comp::device& dev : comp::system::devices())
	{
		_Dev _dev;
		_dev.dev = dev;
		_dev.devtype = dev.type();

		std::string strExtensions = dev.get_info<CL_DEVICE_EXTENSIONS>();

		// needs double type support?
		if(std::is_same<t_real, double>::value)
		{
			bool bHasDouble = (strExtensions.find("cl_khr_fp64") != std::string::npos);
			if(!bHasDouble) continue;
		}

		vecAllDevs.push_back(_dev);
	}

	if(vecAllDevs.size() == 0)
		return boost::optional<comp::device>();


	// assign and sort by device score
	std::sort(vecAllDevs.begin(), vecAllDevs.end(),
		[](const _Dev& dev1, const _Dev& dev2) -> bool
		{
			int (*get_device_score)(cl_device_type ty) = [](cl_device_type ty) -> int
			{
				int iScore = 0;

				if(ty & CL_DEVICE_TYPE_GPU)
					iScore += 1000;
				if(ty & CL_DEVICE_TYPE_ACCELERATOR)
					iScore += 100;
				if(ty & CL_DEVICE_TYPE_CPU)
					iScore += 10;

				return iScore;
			};

			int iScore1 = get_device_score(dev1.devtype);
			int iScore2 = get_device_score(dev2.devtype);

			return iScore1 > iScore2;
		});

	return boost::optional<comp::device>(vecAllDevs[0].dev);
}


/**
 * add CL type definitions to a Boost.Comp function
 */
template<typename t_scalar, class COMP_FKT>
void inject_comp_typedefs(COMP_FKT& fkt)
{
	std::string strSrc = fkt.source();

	// start of function block
	std::size_t iBlockStart = strSrc.find('{');
	if(iBlockStart == std::string::npos)
		return;

	// type def string
	const std::string& strTypeDefs = get_cl_typedefs<t_scalar>();

	strSrc.insert(iBlockStart+1, "\n" + strTypeDefs + "\n");
	fkt.set_source(strSrc);
}


}
#endif
