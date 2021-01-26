/**
 * tlibs test file
 * @author Tobias Weber <tobias.weber@tum.de>
 * @license GPLv2 or GPLv3
 */

// clang -DNO_IOSTR -o loaddat loaddat.cpp ../log/log.cpp -std=c++11 -lstdc++ -lboost_iostreams -lm

#include <iostream>
#include "../file/loaddat.h"

int main()
{
	//tl::comp_stream_to_stream<char>(std::cin, std::cout);
	//tl::comp_stream_to_stream<wchar_t>(std::wcin, std::wcout);

	//std::wstring str(L"123");
	//std::wcout << tl::str_to_var<int, std::wstring>(str) << std::endl;

	//tl::DatFile<float, wchar_t> dat;

	tl::DatFile<float, char> dat;
	//tl::DatFile<float, wchar_t> dat;
	dat.SetSeparatorChars({':'});
	//dat.Load("/home/tw/Measurements/mira-mnsi-14/data/9925_00009990.dat");
	//dat.Load(L"/home/tw/Measurements/mira-mnsi-14/data/9925_00009990.dat");
	dat.Load("/home/tweber/Messdaten/mira-mnsi-15/data/10365_00015382.dat");
	if(dat.IsOk())
	{
		std::cout << "Number of columns: " << dat.GetColumnCount() << std::endl;
		std::cout << "Number of rows: " << dat.GetColumn(0).size() << std::endl;

		const auto& hdr = dat.GetHeader();
		for(const auto& pair : hdr)
		{
			std::cout << pair.first << " = " << pair.second << std::endl;
			//std::wcout << pair.first << L" = " << pair.second << std::endl;
		}

		dat.Save("/tmp/tst.dat");
	}

	return 0;
}
