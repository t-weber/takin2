/**
 * tlibs test file
 * @author Tobias Weber <tobias.weber@tum.de>
 * @license GPLv2 or GPLv3
 */

// clang -o typeidx typeidx.cpp -std=c++11 -lstdc++ -lm

#include "../log/debug.h"
#include <iostream>
#include <vector>

template<class T> void tst(T&& t)
{
	std::cout << tl::get_typename<decltype(t)>() << std::endl;
}

void tst2(int&& i) {}


int main(int argc, char** argv)
{
	std::vector<int> a = {1,2,3,4,5};
	const std::vector<int>& b = a;

	std::vector<bool> c;

	std::cout << tl::get_typename<decltype(main)>() << std::endl;
	std::cout << tl::get_typename<decltype(&main)>() << std::endl;
	std::cout << std::endl;

	std::cout << tl::get_typename<decltype(a)>() << std::endl;
	std::cout << tl::get_typename<decltype(a[0])>() << std::endl;
	std::cout << tl::get_typename<decltype(b)>() << std::endl;
	std::cout << std::endl;

	std::cout << tl::get_typename<decltype(c)>() << std::endl;
	std::cout << tl::get_typename<decltype(c[0])>() << std::endl;
	std::cout << std::endl;


	int i=0;
	tst(i);
	tst(int(i));

	//tst2(i);
	tst2(int(i));

	tst(c);
	tst(std::move(c));
	tst(std::forward<decltype(c)>(c));
	tst(std::forward<decltype(b)>(b));
	std::cout << std::endl;

	return 0;
}
