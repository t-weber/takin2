/**
 * random distributions
 * @author Tobias Weber <tobias.weber@tum.de>
 * @date sep-2016
 * @license GPLv2 or GPLv3
 */

#include "distr.h"
#include "log/log.h"
#include <memory>

int main(int argc, char** argv)
{
	try
	{
		std::unique_ptr<QApplication> pApp = std::make_unique<QApplication>(argc, argv);
		::setlocale(LC_ALL, "C");
		QLocale::setDefault(QLocale::English);

		std::unique_ptr<DistrDlg> pDlg = std::make_unique<DistrDlg>();
		pDlg->show();
		return pApp->exec();
	}
	catch(const std::exception& ex)
	{
		tl::log_err(ex.what());
	}
	return -1;
}
