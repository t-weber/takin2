/**
 * mini scan monitor
 * @author Tobias Weber <tobias.weber@tum.de>
 * @date 2015
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

#include "dialog.h"
#include <boost/iostreams/device/file_descriptor.hpp>
#include <boost/iostreams/stream.hpp>

namespace ios = boost::iostreams;


static FILE *pipeProg = nullptr;
static ios::file_descriptor_sink *pfds = nullptr;
static ios::stream_buffer<ios::file_descriptor_sink> *psbuf = nullptr;
static std::ostream *postrProg = nullptr;


bool open_progress(const std::string& strTitle)
{
	close_progress();

	pipeProg = (FILE*)popen(("dialog --title \"" + strTitle + "\" --gauge \"\" 10 50 0").c_str(), "w");
	if(!pipeProg) return false;
	pfds = new ios::file_descriptor_sink(fileno(pipeProg), ios::close_handle);
	psbuf = new ios::stream_buffer<ios::file_descriptor_sink>(*pfds);
	postrProg = new std::ostream(psbuf);

	return true;
}


void close_progress()
{
	if(pipeProg) { pclose(pipeProg); pipeProg = nullptr; }
	if(postrProg) { delete postrProg; postrProg = nullptr; }
	if(psbuf) { delete psbuf; psbuf = nullptr; }
	if(pfds) { delete pfds; pfds = nullptr; }
}


void set_progress(int iPerc, const std::string& strTxt)
{
	if(!postrProg) return;
	(*postrProg) << "XXX\n" << iPerc << "\n"
		<< strTxt << "\nXXX"
		<< std::endl;
}


/*
// g++ -o dialog dialog_indirect.cpp -std=c++11 -lboost_iostreams

#include <chrono>
#include <thread>

int main(int argc, char** argv)
{
	open_progress("Test");
	set_progress(10, "abc\ndef\nghi");
	std::this_thread::sleep_for(std::chrono::milliseconds(250));
	set_progress(20, "xyz\n123\n456");
	std::this_thread::sleep_for(std::chrono::milliseconds(250));

	close_progress();
	return 0;
}*/
