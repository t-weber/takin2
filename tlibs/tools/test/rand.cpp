/**
 * tlibs test file
 * @author Tobias Weber <tobias.weber@tum.de>
 * @license GPLv2 or GPLv3
 */

// g++ -std=c++11 -o rand rand.cpp ../math/rand.cpp ../log/log.cpp -lm -lpthread

#include "../math/rand.h"
#include "../math/stat.h"
#include "../string/string.h"

#include <iostream>
#include <iomanip>
#include <thread>
#include <chrono>
#include <boost/type_traits/function_traits.hpp>


using t_real = double;
using t_seed = typename boost::function_traits<decltype(tl::init_rand_seed)>::arg1_type;


std::mutex mtx;

void rnd_fkt(t_seed iSeed)
{
	mtx.lock();
	tl::init_rand_seed(iSeed);

	std::cout << "seed: " << iSeed << std::endl;
	std::cout << tl::rand01<float>() << " " << tl::rand01<float>() << std::endl;
	std::cout << tl::rand_minmax<float>(0., 10.) << std::endl;
	std::cout << tl::rand_minmax<int>(0, 10) << std::endl;
	std::cout << tl::rand_binomial<int, float>(100, 0.5) << std::endl;

	auto vecRnd = tl::rand_norm_nd<>({1., 2., 3.}, {0.25, 0.5, 0.75});
	for(auto d : vecRnd)
		std::cout << d << ", ";
	std::cout << std::endl;

	vecRnd = tl::rand_norm_nd<>({1., 2., 3.}, {0.25, 0.5, 0.75});
	for(auto d : vecRnd)
		std::cout << d << ", ";
	std::cout << std::endl;

	auto vecRnd2 = tl::rand_exp_nd<>({1., 2., 3.});
	for(auto d : vecRnd2)
		std::cout << d << ", ";
	std::cout << std::endl;

	mtx.unlock();
}




std::tuple<std::vector<t_real>, std::vector<t_real>, std::vector<t_real>, std::vector<t_real>, std::vector<t_real>>
rnd2(std::size_t N, t_real p=0.5, t_real l1=1., t_real l2=-1.)
{
	t_real x{0};
	std::vector<t_real> ns, xs, xs_mean, xs_std, xs2_mean;
	ns.reserve(N); xs.reserve(N); xs_mean.reserve(N); xs_std.reserve(N); xs2_mean.reserve(N);

	for(std::size_t n=0; n<N; ++n)
	{
		t_real r01 = tl::rand01<t_real>();
		x += (r01 < p ? l1 : l2);

		ns.push_back(n);
		xs.push_back(x);
		xs_mean.push_back(tl::mean_value(xs));			// <x>
		xs2_mean.push_back(tl::mean_square_value(xs));	// <x^2>
		//xs_std.push_back(tl::std_dev(xs, 0));			// sqrt(<x^2> - <x>^2)
		xs_std.push_back(std::sqrt(*xs2_mean.rbegin() - *xs_mean.rbegin()* *xs_mean.rbegin()));		// sqrt(<x^2> - <x>^2)
	}

	return std::make_tuple(ns, xs, xs_mean, xs_std, xs2_mean);
}




std::tuple<std::vector<t_real>, std::vector<t_real>>
rnd3(std::size_t N, t_real inc=1., t_real max=1.)
{
	std::vector<t_real> Ns;
	std::vector<t_real> Ks;
	std::vector<t_real> means;

	Ks.reserve(N);
	means.reserve(N);

	for(std::size_t n=0; n<N; ++n)
	{
		t_real k{0};
		t_real x{0};

		while(1)
		{
			x += tl::rand01<t_real>();
			k += inc;

			if(x > max)
				break;
		}

		Ns.push_back(n);
		Ks.push_back(k);
		means.push_back(tl::mean_value(Ks));
	}

	return std::make_tuple(Ns, means);
}



int main(int argc, char** argv)
{
	// test 1
	{
		std::cout << "--------------------------------------------------------------------------------" << std::endl;
		std::cout << "Test 1" << std::endl;
		std::cout << "--------------------------------------------------------------------------------" << std::endl;

		std::thread th1([]{ rnd_fkt(1); });
		std::thread th2([]{ rnd_fkt(2); });

		th1.join();
		th2.join();

		std::cout << "--------------------------------------------------------------------------------\n" << std::endl;
	}



	std::size_t N = 64;
	if(argc > 1)
		N = tl::str_to_var<std::size_t>(std::string{argv[1]});


	// init seed to seconds since epoch
	tl::init_rand_seed(std::chrono::duration_cast<std::chrono::duration<t_seed, std::ratio<1, 1>>>(
		std::chrono::system_clock::now().time_since_epoch()).count());



	// test 2
	{
		std::cout << "--------------------------------------------------------------------------------" << std::endl;
		std::cout << "Test 2" << std::endl;
		std::cout << "--------------------------------------------------------------------------------" << std::endl;

		std::size_t skip = 0;
		if(argc > 2)
			skip = tl::str_to_var<std::size_t>(std::string{argv[2]});

		std::vector<t_real> ns, xs, xs_mean, xs_std, xs2_mean;
		std::tie(ns, xs, xs_mean, xs_std, xs2_mean) = rnd2(N, 0.5, 1., -1.);

		std::cout.precision(5);
		std::cout << "# "
			<< std::setw(10) << "n" << " "
			<< std::setw(10) << "x" << " "
			<< std::setw(10) << "<x>" << " "
			<< std::setw(10) << "std_dev" << " "
			<< std::setw(10) << "<x^2>" << "\n";

		for(std::size_t i=0; i<ns.size(); ++i)
		{
			if(skip && i%skip!=0)
				continue;

			std::cout << "  "
				<< std::setw(10) << ns[i] << " "
				<< std::setw(10) << xs[i] << " "
				<< std::setw(10) << xs_mean[i] << " "
				<< std::setw(10) << xs_std[i] << " "
				<< std::setw(10) << xs2_mean[i] << "\n";
		}

		std::cout << "--------------------------------------------------------------------------------\n" << std::endl;
	}



	// test 3
	{
		std::cout << "--------------------------------------------------------------------------------" << std::endl;
		std::cout << "Test 3" << std::endl;
		std::cout << "--------------------------------------------------------------------------------" << std::endl;

		std::vector<t_real> ns, means;
		std::tie(ns, means) = rnd3(N, 1., 1.);

		std::cout.precision(5);
		std::cout << "# "
			<< std::setw(10) << "n" << " "
			<< std::setw(10) << "<k>" << "\n";

		for(std::size_t i=0; i<ns.size(); ++i)
		{
			std::cout << "  "
				<< std::setw(10) << ns[i] << " "
				<< std::setw(10) << means[i] << "\n";
		}

		std::cout << "--------------------------------------------------------------------------------\n" << std::endl;
	}

	return 0;
}
