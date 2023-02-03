/**
 * Parse & Runtime Info
 * @author Tobias Weber <tobias.weber@tum.de>
 * @date 12-oct-14
 * @license GPLv2 or GPLv3
 *
 * ----------------------------------------------------------------------------
 * tlibs -- a physical-mathematical C++ template library
 * Copyright (C) 2017-2021  Tobias WEBER (Institut Laue-Langevin (ILL),
 *                          Grenoble, France).
 * Copyright (C) 2015-2017  Tobias WEBER (Technische Universitaet Muenchen
 *                          (TUM), Garching, Germany).
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 2 of the License, or
 * (at your option) version 3.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 * ----------------------------------------------------------------------------
 */

#include "info.h"

ParseInfo::ParseInfo()
{
	pmapModules = new t_mods();
	pGlobalSyms = new SymbolTable();
	phandles = new HandleManager();
	pmutexGlobal = new std::mutex();
	pmutexGlobalSyms = new std::mutex();
	pmutexTraceback = new std::mutex();
}

ParseInfo::~ParseInfo()
{
	//if(phandles) { delete phandles; phandles=0; }
	if(pGlobalSyms) { delete pGlobalSyms; pGlobalSyms=0; }
	if(pmutexGlobal) { delete pmutexGlobal; pmutexGlobal=0; }
	if(pmutexGlobalSyms) { delete pmutexGlobalSyms; pmutexGlobalSyms=0; }
	if(pmutexTraceback) { delete pmutexTraceback; pmutexTraceback=0; }

	if(pmapModules)
	{
		for(ParseInfo::t_mods::value_type vals : *pmapModules)
		{
			if(vals.second)
				delete vals.second;
		}

		pmapModules->clear();
		delete pmapModules;
		pmapModules = 0;
	}
}

void ParseInfo::PushTraceback(std::string&& strTrace)
{
	if(!bEnableDebug) return;

	t_oneTraceback *pStck = 0;
	{
		std::lock_guard<std::mutex> lck(*pmutexTraceback);
		pStck = &stckTraceback[std::this_thread::get_id()];
	}

	pStck->push_front(std::forward<std::string>(strTrace));
}

void ParseInfo::PopTraceback()
{
	if(!bEnableDebug) return;

	t_oneTraceback *pStck = 0;
	{
		std::lock_guard<std::mutex> lck(*pmutexTraceback);
		pStck = &stckTraceback[std::this_thread::get_id()];
	}

	pStck->pop_front();
}
