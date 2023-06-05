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
#include <dialog/dialog.h>


static void *g_pProgress = nullptr;
std::string g_strTitle;


bool open_progress(const std::string& strTitle)
{
	g_strTitle = strTitle;
	close_progress();

	init_dialog(stdin, stdout);
	dlg_clear();

	g_pProgress = dlg_allocate_gauge(g_strTitle.c_str(), "", 10, 50, 0);
	if(g_pProgress)
		return true;
	return false;
}


void close_progress()
{
	if(g_pProgress)
	{
		dlg_free_gauge(g_pProgress);
		g_pProgress = nullptr;
	}
	end_dialog();
}


void set_progress(int iPerc, const std::string& strTxt)
{
	if(!g_pProgress) return;

	g_pProgress = dlg_reallocate_gauge(g_pProgress, g_strTitle.c_str(), strTxt.c_str(), 10, 50, iPerc);
	dlg_update_gauge(g_pProgress, iPerc);
}


/*
// g++ -o dialog dialog.cpp -std=c++11 -ldialog

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
