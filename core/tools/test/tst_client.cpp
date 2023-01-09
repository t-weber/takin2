/**
 * test instrument client
 * @author Tobias Weber <tobias.weber@tum.de>
 * @date 2014
 * @license GPLv2
 * clang -o tst_server -I. -I../.. ../../tools/misc/tst_server.cpp ../../tlibs/net/tcp.cpp ../../tlibs/log/log.cpp -lstdc++ -std=c++11 -lboost_system -lboost_iostreams -lpthread -lm
 *
 * ----------------------------------------------------------------------------
 * Takin (inelastic neutron scattering software package)
 * Copyright (C) 2017-2021  Tobias WEBER (Institut Laue-Langevin (ILL),
 *                          Grenoble, France).
 * Copyright (C) 2013-2017  Tobias WEBER (Technische Universitaet Muenchen
 *                          (TUM), Garching, Germany).
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; version 2 of the License.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along
 * with this program; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 * ----------------------------------------------------------------------------
 */

/*
 * clang -o tst_client -I../.. tst_client.cpp ../../tlibs/net/tcp.cpp ../../tlibs/log/log.cpp -lstdc++ -std=c++11 -lboost_system -lpthread -lm
 */

#include <fstream>
#include "../../tlibs/net/tcp.h"
#include "../../tlibs/log/log.h"

using namespace tl;

struct TstOut
{
	std::ofstream ofstr;

	TstOut() : ofstr("tst_out.txt")
	{}

	void print(const std::string& str)
	{
		//ofstr << str << std::endl;
		std::cout << "out: " << str << std::endl;
	}
};

static void disconnected(const std::string& strHost, const std::string& strSrv)
{
	log_info("Disconnected from ", strHost, " on port ", strSrv, ".");
}

static void connected(const std::string& strHost, const std::string& strSrv)
{
	log_info("Connected to ", strHost, " on port ", strSrv, ".");
}

static void received(const std::string& strMsg)
{
	log_info("Received: ", strMsg, ".");
}


int main(int argc, char** argv)
{
	if(argc < 3)
	{
		std::cerr << "Usage: " << argv[0] << " <server> <port>" << std::endl;
		return -1;
	}


	TcpTxtClient<> client;
	//TstOut tstout;
	//client.add_receiver(boost::bind(&TstOut::print, &tstout, _1));
	client.add_receiver(received);
	client.add_disconnect(disconnected);
	client.add_connect(connected);


	if(!client.connect(argv[1], argv[2]))
	{
		log_err("Cannot connect.");
		return -1;
	}

	std::string strMsg;
	while(client.is_connected())
	{
		std::getline(std::cin, strMsg);
		if(strMsg == "!exit!")
			break;

		strMsg+="\n";
		client.write(strMsg);
	}

	return 0;
}
