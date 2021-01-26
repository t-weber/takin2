/**
 * invocation of gnuplot
 * @autor Tobias Weber <tobias.weber@tum.de>
 * @date 24-dec-2013
 * @license GPLv2 or GPLv3
 */

#include "gnuplot.h"
#include "gnuplot_impl.h"

namespace tl
{
	template struct PlotObj<double>;
	template struct PlotObj<float>;

	template class GnuPlot<double>;
	template class GnuPlot<float>;
}
