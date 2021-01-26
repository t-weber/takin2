/**
 * tlibs test file
 * @author Tobias Weber <tobias.weber@tum.de>
 * @license GPLv2 or GPLv3
 */

// gcc -o dft3 dft3.cpp ../math/fftw.cpp ../log/log.cpp -lfftw3 -lstdc++ -lm -std=c++11

#include "../math/dft.h"
#include "../math/fftw.h"
#include "../time/stopwatch.h"
#include <iostream>

#define NUM 4096

int main()
{
	double eps = 1e-4;

	double* dInR = new double[NUM];
	double* dInI = new double[NUM];

	for(unsigned int i=0; i<NUM; ++i)
	{
		dInR[i] = i+1;
		dInI[i] = NUM-i;
	}

	const std::size_t iSize = NUM;

	double *dOutR1 = new double[iSize], *dOutI1 = new double[iSize];
	double *dOutR2 = new double[iSize], *dOutI2 = new double[iSize];
	double *dOutR3 = new double[iSize], *dOutI3 = new double[iSize];

	tl::FFT<double> fourier1(iSize);
	tl::DFT<double> fourier2(iSize);
	tl::FFTw fourier3(iSize);

	tl::Stopwatch<double> sw1, sw2, sw3;

	std::cout << "Forward trafo 1..." << std::endl;
	sw1.start(); fourier1.trafo(dInR, dInI, dOutR1, dOutI1, 0); sw1.stop();
	std::cout << "Forward trafo 2..." << std::endl;
	sw2.start(); fourier2.trafo(dInR, dInI, dOutR2, dOutI2, 0); sw2.stop();
	std::cout << "Forward trafo 3..." << std::endl;
	sw3.start(); fourier3.trafo(dInR, dInI, dOutR3, dOutI3, 0); sw3.stop();

	for(int i=0; i<iSize; ++i)
	{
		if(!tl::float_equal(dOutR1[i], dOutR2[i], eps)
			|| !tl::float_equal(dOutR1[i], dOutR3[i], eps)
			|| !tl::float_equal(dOutI1[i], dOutI2[i], eps)
			|| !tl::float_equal(dOutI1[i], dOutI3[i], eps))
		{
			std::cerr << "Mismatch at i=" << i << "!" << std::endl;
			std::cout << "1: " << dOutR1[i] << " + i*" << dOutI1[i] << std::endl;
			std::cout << "2: " << dOutR2[i] << " + i*" << dOutI2[i] << std::endl;
			std::cout << "3: " << dOutR3[i] << " + i*" << dOutI3[i] << std::endl;
			return -1;
		}
	}

	std::cout << "Times needed: \n"
		<< "\tfft: " << sw1.GetDur() << "s\n"
		<< "\tdft: " << sw2.GetDur() << "s\n"
		<< "\tfftw: " << sw3.GetDur() << "s\n";

	std::cout << std::endl;

/*	std::cout << "fft: ";
	for(int i=0; i<iSize; ++i)
	{
		std::cout << dOutR1[i];
		if(!tl::float_equal(dOutI1[i], 0., 1e-7))
			std::cout << " + " << dOutI1[i] << "*i";
		std::cout << ", ";
	}
	std::cout << std::endl;

	std::cout << "dft: ";
	for(int i=0; i<iSize; ++i)
	{
		std::cout << dOutR2[i];
		if(!tl::float_equal(dOutI2[i], 0., 1e-7))
			std::cout << " + " << dOutI2[i] << "*i";
		std::cout << ", ";
	}
	std::cout << std::endl;

	std::cout << "fftw: ";
	for(int i=0; i<iSize; ++i)
	{
		std::cout << dOutR3[i];
		if(!tl::float_equal(dOutI3[i], 0., 1e-7))
			std::cout << " + " << dOutI3[i] << "*i";
		std::cout << ", ";
	}
	std::cout << std::endl;*/



	//std::vector<std::complex<double>> vec = tl::arrs_to_cvec(dOutR, dOutI, iSize);
	//vec = tl::dft_shift(vec, 1.);
	//tl::cvec_to_arrs(vec, dOutR, dOutI);




	double *dOutR4 = new double[iSize], *dOutI4 = new double[iSize];
	double *dOutR5 = new double[iSize], *dOutI5 = new double[iSize];
	double *dOutR6 = new double[iSize], *dOutI6 = new double[iSize];

	std::cout << "Backward trafo 1..." << std::endl;
	sw1.start(); fourier1.trafo(dOutR1, dOutI1, dOutR4, dOutI4, 1); sw1.stop();
	std::cout << "Backward trafo 2..." << std::endl;
	sw2.start(); fourier2.trafo(dOutR2, dOutI2, dOutR5, dOutI5, 1); sw2.stop();
	std::cout << "Backward trafo 3..." << std::endl;
	sw3.start(); fourier3.trafo(dOutR3, dOutI3, dOutR6, dOutI6, 1); sw3.stop();



	for(int i=0; i<iSize; ++i)
	{
		if(!tl::float_equal(dOutR4[i], dOutR5[i], eps)
			|| !tl::float_equal(dOutR4[i], dOutR6[i], eps)
			|| !tl::float_equal(dOutI4[i], dOutI5[i], eps)
			|| !tl::float_equal(dOutI4[i], dOutI6[i], eps))
		{
			std::cerr << "Mismatch at i=" << i << "!" << std::endl;
			std::cout << "4: " << dOutR4[i] << " + i*" << dOutI4[i] << std::endl;
			std::cout << "5: " << dOutR5[i] << " + i*" << dOutI5[i] << std::endl;
			std::cout << "6: " << dOutR6[i] << " + i*" << dOutI6[i] << std::endl;
			return -1;
		}
	}


	std::cout << "Times needed: \n"
		<< "\tfft: " << sw1.GetDur() << "s\n"
		<< "\tdft: " << sw2.GetDur() << "s\n"
		<< "\tfftw: " << sw3.GetDur() << "s\n";

	std::cout << std::endl;


/*	std::cout << "ifft: ";
	for(int i=0; i<iSize; ++i)
	{
		std::cout << dOutR4[i]/iSize;
		if(!tl::float_equal(dOutI4[i], 0., 1e-7))
			std::cout << " + " << dOutI4[i]/iSize << "*i";
		std::cout << ", ";
	}
	std::cout << std::endl;

	std::cout << "idft: ";
	for(int i=0; i<iSize; ++i)
	{
		std::cout << dOutR5[i]/iSize;
		if(!tl::float_equal(dOutI5[i], 0., 1e-7))
			std::cout << " + " << dOutI5[i]/iSize << "*i";
		std::cout << ", ";
	}
	std::cout << std::endl;

	std::cout << "ifftw: ";
	for(int i=0; i<iSize; ++i)
	{
		std::cout << dOutR6[i]/iSize;
		if(!tl::float_equal(dOutI6[i], 0., 1e-7))
			std::cout << " + " << dOutI6[i]/iSize << "*i";
		std::cout << ", ";
	}
	std::cout << std::endl;*/


	delete[] dInR; delete[] dInI;
	delete[] dOutR1; delete[] dOutI1;
	delete[] dOutR2; delete[] dOutI2;
	delete[] dOutR3; delete[] dOutI3;
	delete[] dOutR4; delete[] dOutI4;
	delete[] dOutR5; delete[] dOutI5;
	delete[] dOutR6; delete[] dOutI6;

	return 0;
}
