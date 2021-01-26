/**
 * tlibs test file
 * Test distributions
 * @author Tobias Weber <tobias.weber@tum.de>
 * @date 10-sep-2016
 * @license GPLv2 or GPLv3
 */

// gcc -std=c++11 -o distr distr.cpp ../gfx/gnuplot.cpp ../log/log.cpp -lboost_iostreams -lstdc++ -lm

#include <iostream>
#include <vector>
#include <boost/math/distributions/binomial.hpp>
#include <boost/math/distributions/hypergeometric.hpp>
#include <boost/math/distributions/poisson.hpp>
#include <boost/math/distributions/normal.hpp>
#include <boost/math/distributions/students_t.hpp>
#include <boost/math/distributions/chi_squared.hpp>
#include <boost/math/distributions/cauchy.hpp>
#include "../gfx/gnuplot.h"
#include "../math/math.h"


namespace m = boost::math;

using t_real = double;
//using t_dist = m::normal_distribution<t_real>;
//using t_dist = m::students_t_distribution<t_real>;
//using t_dist = m::binomial_distribution<t_real>;
//using t_dist = m::hypergeometric_distribution<t_real>;
//using t_dist = m::poisson_distribution<t_real>;
//using t_dist = m::chi_squared_distribution<t_real>;
using t_dist = m::cauchy_distribution<t_real>;
//#define SHOW_PROPERTIES


int main()
{
	//t_dist dist(2.5, 1.);		// normal
	//t_dist dist(1.);			// t
	//t_dist dist(5., 0.5);		// binomial
	//t_dist dist(2., 5., 10.);	// hypergeo
	//t_dist dist(3.);			// poisson
	//t_dist dist(4.);			// chi^2
	t_dist dist(2.5, 0.5);		// Cauchy / Lorentzian


#ifdef SHOW_PROPERTIES
	try
	{
		std::cout
			<< "mean: " << m::mean(dist)
			<< ", median: " << m::median(dist)
			<< ", mode: " << m::mode(dist)
			<< ", kurtosis: " << m::kurtosis(dist)
			<< ", kurtosis ex: " << m::kurtosis_excess(dist)
			<< ", skewness: " << m::skewness(dist)
			<< ", std dev: " << m::standard_deviation(dist)
			<< std::endl;

		std::cout << "area (only if density function is symmetric): " << m::cdf(dist, m::mean(dist))*2.
			<< std::endl;
	}
	catch(const std::exception& ex)
	{
		std::cerr << ex.what() << std::endl;
	}
#endif

	tl::PlotObj<t_real> line1, line2, line3;
	line1.linestyle = tl::STYLE_LINES_SOLID;
	line2.linestyle = tl::STYLE_LINES_SOLID;
	line3.linestyle = tl::STYLE_LINES_SOLID;
	line1.odSize = 2.;
	line2.odSize = 2.;
	line3.odSize = 2.;
	line1.oiColor = 0xff0000;
	line2.oiColor = 0x0000ff;
	line3.oiColor = 0x007700;

	line1.vecX = line2.vecX = tl::linspace(0., 5., 128);
	line3.vecX = tl::linspace(0.1, 0.9, 128);
	for(t_real x : line1.vecX)
	{
		line1.vecY.push_back(m::cdf(dist, x));
		line2.vecY.push_back(m::pdf(dist, x));
	}
	for(t_real x : line3.vecX)
	{
		line3.vecY.push_back(m::quantile(dist, x));
	}

	tl::GnuPlot<t_real> gpl;
	gpl.Init();

	gpl.SetTerminal(0, "wxt");
	gpl.StartPlot();
	gpl.AddLine(line1);
	gpl.AddLine(line2);
	gpl.FinishPlot();

	gpl.SetTerminal(1, "wxt");
	gpl.StartPlot();
	gpl.AddLine(line3);
	gpl.FinishPlot();

	return 0;
}
