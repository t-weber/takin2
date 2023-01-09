/**
 * @author Tobias Weber <tobias.weber@tum.de>
 * @license GPLv2
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

// gcc -o tst_proc tst_proc.cpp ../../tlibs/math/rand.cpp ../../tlibs/log/log.cpp -lpthread -lrt -lstdc++ -std=c++11

#include <iostream>
#include <unistd.h>
#include <string>
#include <boost/interprocess/ipc/message_queue.hpp>
#include "../../tlibs/math/rand.h"

namespace proc = boost::interprocess;

enum class ProcMsgTypes
{
	PRINT,
	QUIT
};

struct ProcMsg
{
	ProcMsgTypes ty;

	int iMsg;
};


void child(proc::message_queue& msgToParent, proc::message_queue& msgFromParent)
{
	std::cout << "In child process." << std::endl;

	while(1)
	{
		ProcMsg themsg;
		std::size_t iSize;
		unsigned int iPrio;
		msgFromParent.receive(&themsg, sizeof(themsg), iSize, iPrio);

		if(iSize != sizeof(themsg))
			std::cerr << "Size mismatch." << std::endl;

		switch(themsg.ty)
		{
			case ProcMsgTypes::PRINT:
			{
				std::cout << "Child received: " << themsg.iMsg << std::endl;

				ProcMsg themsgOut;
				themsgOut.iMsg = themsg.iMsg * themsg.iMsg;
				msgToParent.send(&themsgOut, sizeof(themsgOut), 0);
				break;
			}
			case ProcMsgTypes::QUIT:
			{
				std::cout << "Exiting child process." << std::endl;
				exit(0);
				break;
			}
		}
	}
}


void parent(proc::message_queue& msgFromChild, proc::message_queue& msgToChild)
{
	std::cout << "In parent process." << std::endl;

	ProcMsg themsgIn;
	std::size_t iSizeIn;
	unsigned int iPrioIn;

	ProcMsg themsg;
	themsg.ty = ProcMsgTypes::PRINT;

	for(int i=0; i<10; ++i)
	{
		themsg.iMsg = i+1;
		msgToChild.send(&themsg, sizeof(themsg), 0);

		msgFromChild.receive(&themsgIn, sizeof(themsgIn), iSizeIn, iPrioIn);
		std::cout << "Parent received: " << themsgIn.iMsg << std::endl;
	}

	themsg.ty = ProcMsgTypes::QUIT;
	msgToChild.send(&themsg, sizeof(themsg), 0);
}


int main()
{
	try
	{
		tl::init_rand();
		std::string strName = std::string("tst_proc_") + tl::rand_name<std::string>(8);
		std::cout << "Process name: " << strName << std::endl;

		proc::message_queue msgOut(proc::create_only, (strName + "_out").c_str(), 16, sizeof(ProcMsg));
		proc::message_queue msgIn(proc::create_only, (strName + "_in").c_str(), 16, sizeof(ProcMsg));

		pid_t pid = fork();
		if(pid == 0)
		{
			child(msgIn, msgOut);
		}
		else
		{
			std::cout << "Forked child: " << pid << "." << std::endl;
			parent(msgIn, msgOut);

			std::cout << "Exiting parent process." << std::endl;
			proc::message_queue::remove((strName + "_out").c_str());
			proc::message_queue::remove((strName + "_in").c_str());
		}
	}
	catch(const std::exception& ex)
	{
		std::cerr << ex.what() << std::endl;
	}

	return 0;
}
