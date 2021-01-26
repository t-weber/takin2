/**
 * tlibs test file
 * @author Tobias Weber <tobias.weber@tum.de>
 * @license GPLv2 or GPLv3
 */

// clang -o thread thread.cpp -std=c++11 -lstdc++ -lpthread

#include "../helper/thread.h"
#include <iostream>

int tst(int i)
{
	return i*i;
}

int main()
{
	tl::ThreadPool<int()> tp;
	for(int i=0; i<10; ++i)
		tp.AddTask(std::bind(tst, i));

	for(int i=0; i<10; ++i)
		tp.AddTask([i]()->int{ return i; });

	tp.StartTasks();

	auto& lstFut = tp.GetFutures();
	for(auto& fut : lstFut)
		std::cout << fut.get() << std::endl;
	return 0;
}
