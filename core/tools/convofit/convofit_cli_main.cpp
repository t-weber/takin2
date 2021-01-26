/**
 * Convolution fitting
 * @author Tobias Weber <tobias.weber@tum.de>
 * @date dec-2015
 * @license GPLv2
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
