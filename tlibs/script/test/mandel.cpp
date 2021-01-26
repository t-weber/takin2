// clang -march=native -O3 -o mandel mandel.cpp -std=c++11 -lstdc++ -lm
// hermelin performance comparison

#include <iostream>
#include <fstream>
#include <iomanip>
#include <complex>
#include <cstdlib>

typedef std::complex<double> t_cplx;

unsigned int mandel(const t_cplx& c)
{
	t_cplx z(0., 0.);

	unsigned int iIter=0;
	for(iIter=0; iIter<128; ++iIter)
	{
		if(std::norm(z) > 4.)
			break;

		z = z*z + c;
	}

	return iIter;
}


int main()
{
	const int iXSize = 512;
	const int iYSize = 512;

	double dXScale = 1./iXSize * 2.;
	double dYScale = 1./iYSize * 2.;
	double dXOffs = -0.5;
	double dYOffs = 0.;


	std::ofstream ofstr("mandel.dat");

	for(int iY=0; iY<iYSize; ++iY)
	{
		for(int iX=0; iX<iXSize; ++iX)
		{
			t_cplx c((iX - iXSize/2)*dXScale + dXOffs,
				(iY - iYSize/2)*dYScale + dYOffs);
			unsigned int iIter = mandel(c);

			ofstr << std::setw(16) << iIter;
		}
		ofstr << "\n";

		std::cout << "Line " << iY << std::endl;
	}

	ofstr.flush();
	ofstr.close();


	std::system(("gnuplot -p -e \"set xrange [0:" + std::to_string(iXSize) + "];" +
			"set yrange [0:" + std::to_string(iYSize) + "];" +
			"plot 'mandel.dat' matrix with image\n\"").c_str());
	return 0;
}
