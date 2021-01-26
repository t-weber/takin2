/**
 * tlibs test file
 * @author Tobias Weber <tobias.weber@tum.de>
 * @license GPLv2 or GPLv3
 */

// clang -o ferro ferro.cpp ../gfx/gnuplot.cpp ../log/log.cpp -lstdc++ -lm -lboost_iostreams -std=c++11

#include "../math/mag.h"
#include "../math/neutrons.h"
#include "../gfx/gnuplot.h"

#include <tuple>
#include <iostream>

using namespace tl;
template<class T> using t_magatompos = std::tuple<ublas::vector<T>, std::complex<T>>;

int main()
{
	double k = co::k_B / one_meV * kelvin;
	double a = 1.;
	double J1 = 0.5*k;
	double J2 = 0.5*k;
	double S = 7./2.;

	std::vector<ublas::vector<double>> vecAtomsN0 = get_neighbour_atoms(UCType::FCC, 0, a);
	std::vector<ublas::vector<double>> vecAtomsN1 = get_neighbour_atoms(UCType::FCC, 1, a);

	std::vector<t_magatompos<double>> vecMagAtoms;
	for(const ublas::vector<double>& vec : vecAtomsN0)
	{
		vecMagAtoms.push_back(t_magatompos<double>(vec, J1));
		std::cout << vec << ": " << J1 << std::endl;
	}
	for(const ublas::vector<double>& vec : vecAtomsN1)
	{
		vecMagAtoms.push_back(t_magatompos<double>(vec, J2));
		std::cout << vec << ": " << J2 << std::endl;
	}



	std::vector<double> vecq = linspace(0., 1.5, 128);
	std::vector<double> vecE1, vecE11, vecE111;
	vecE1.reserve(vecq.size());
	vecE11.reserve(vecq.size());
	vecE111.reserve(vecq.size());

	for(double q : vecq)
	{
		ublas::vector<double> vec1 = make_vec({1., 0., 0.});
		ublas::vector<double> vec11 = make_vec({1., 1., 0.}); vec11 /= ublas::norm_2(vec11);
		ublas::vector<double> vec111 = make_vec({1., 1., 1.}); vec111 /= ublas::norm_2(vec111);

		vec1 *= q * 2*M_PI/a; vec11 *= q * 2*M_PI/a; vec111 *= q * 2*M_PI/a;

		double E1 = ferromag(vecMagAtoms, vec1, S);
		double E11 = ferromag(vecMagAtoms, vec11, S);
		double E111 = ferromag(vecMagAtoms, vec111, S);

		vecE1.push_back(E1);
		vecE11.push_back(E11);
		vecE111.push_back(E111);
	}


	GnuPlot<> plt;
	plt.Init();
	plt.SetLegendPlace("bottom right");
	plt.SetXLabel("q (rlu)");
	plt.SetYLabel("E (meV)");
	plt.StartPlot();

	PlotObj<> obj1;
	obj1.linestyle = STYLE_LINES_SOLID;
	obj1.vecX = vecq;
	obj1.vecY = vecE1;
	obj1.oiColor = 0xff0000;
	obj1.strLegend = "[100]";
	plt.AddLine(obj1);

	PlotObj<> obj11;
	obj11.linestyle = STYLE_LINES_SOLID;
	obj11.vecX = vecq;
	obj11.vecY = vecE11;
	obj11.oiColor = 0x009900;
	obj11.strLegend = "[110]";
	plt.AddLine(obj11);

	PlotObj<> obj111;
	obj111.linestyle = STYLE_LINES_SOLID;
	obj111.vecX = vecq;
	obj111.vecY = vecE111;
	obj111.oiColor = 0x0000ff;
	obj111.strLegend = "[111]";
	plt.AddLine(obj111);

	plt.FinishPlot();

	//plt.SimplePlot(vecq, vecE111, std::vector<double>(), std::vector<double>(), STYLE_LINES_SOLID);
	return 0;
}
