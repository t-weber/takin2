/**
 * Convolution fitting
 * @author Tobias Weber <tobias.weber@tum.de>
 * @date dec-2015
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

#include <boost/asio/io_service.hpp>
#include <boost/asio/signal_set.hpp>
#include <boost/scope_exit.hpp>

#include <iostream>
#include <locale>
#include <thread>

#include "convofit_cli.h"
#include "tlibs/log/log.h"

namespace asio = boost::asio;
namespace sys = boost::system;


int main(int argc, char** argv)
{
	std::ios_base::sync_with_stdio(0);

#ifdef NO_TERM_CMDS
	tl::Log::SetUseTermCmds(0);
#endif

	// plain C locale
	/*std::*/setlocale(LC_ALL, "C");
	std::locale::global(std::locale::classic());


	// --------------------------------------------------------------------
	// install exit signal handlers
	asio::io_service ioSrv;
	asio::signal_set sigInt(ioSrv, SIGABRT, SIGTERM, SIGINT);
	sigInt.async_wait([&ioSrv](const sys::error_code& err, int iSig)
	{
		tl::log_warn("Hard exit requested via signal ", iSig, ". This may cause a fault.");
		if(err) tl::log_err("Error: ", err.message(), ", error category: ", err.category().name(), ".");
		ioSrv.stop();
#ifdef SIGKILL
		// TODO: use specific PIDs
		//std::system("killall -s KILL gnuplot");
		std::raise(SIGKILL);
#endif
		exit(-1);
	});
	std::thread thSig([&ioSrv]() { ioSrv.run(); });
	BOOST_SCOPE_EXIT(&ioSrv, &thSig)
	{
		//tl::log_debug("Exiting...");
		ioSrv.stop();
		thSig.join();
	}
	BOOST_SCOPE_EXIT_END
	// --------------------------------------------------------------------


	return convofit_main(argc, argv);
}
