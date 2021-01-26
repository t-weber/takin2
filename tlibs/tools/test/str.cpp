/**
 * tlibs test file
 * @author Tobias Weber <tobias.weber@tum.de>
 * @license GPLv2 or GPLv3
 */

#include "../string/string.h"

int main()
{
	std::cout << tl::var_to_str<int>(1234567890, 10, 3) << std::endl;

	std::vector<int> vec({1,2,3,4,5,6,7,8});
	std::cout << tl::cont_to_str<decltype(vec)>(vec);
	std::cout << std::endl;

	std::cout << tl::str_is_digits(std::string("123456")) << std::endl;
	std::cout << tl::str_is_digits(std::wstring(L"123456")) << std::endl;
	std::cout << tl::str_is_digits(std::string("123a456")) << std::endl;
	std::cout << std::endl;

	std::cout << tl::ends_with<std::string>("Test123", "123") << std::endl;;
	std::cout << tl::ends_with<std::string>("Test1234", "123") << std::endl;;
	std::cout << tl::ends_with<std::string>("Test123", "21") << std::endl;;

	return 0;
}
