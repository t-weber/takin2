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

#ifndef __SCRIPT_INFOS_H__
#define __SCRIPT_INFOS_H__

#include "node.h"
#include "handles.h"
#include "symbol.h"
#include "types.h"
#include "lexer.h"


// stuff that can change during execution
struct RuntimeInfo
{
	// function to execute, e.g. "main" (with external local symbol table)
	t_string strExecFkt;
	t_string strInitScrFile;
	SymbolTable *pLocalSymsOverride = nullptr;

	// currently active function
	const NodeFunction *pCurFunction = nullptr;
	const NodeCall *pCurCaller = nullptr;
	bool bWantReturn = 0;

	const Node* pCurLoop = nullptr;
	bool bWantBreak = 0;
	bool bWantContinue = 0;


	// implicitely return last symbol in function
	bool bImplicitRet = 0;


	bool IsExecDisabled() const
	{
		return bWantReturn || bWantBreak || bWantContinue;
	}
	void EnableExec()
	{
		bWantReturn = bWantBreak = bWantContinue = 0;
	}


	RuntimeInfo() = default;
	RuntimeInfo(const RuntimeInfo&) = delete;
	~RuntimeInfo() = default;
};


// stuff that is (more or less) fixed after parsing
struct ParseInfo
{
	// external imported modules
	typedef std::unordered_map<t_string, Node*> t_mods;
	t_mods *pmapModules = nullptr;

	// all functions from all modules
	typedef std::vector<NodeFunction*> t_funcs;
	t_funcs vecFuncs;

	// global symbol table
	SymbolTable *pGlobalSyms = nullptr;
	std::mutex *pmutexGlobalSyms = nullptr;

	HandleManager *phandles = nullptr;


	// mutex for script if no explicit mutex given
	std::mutex *pmutexGlobal = nullptr;


	bool bEnableDebug = 0;
	std::mutex *pmutexTraceback = nullptr;
	typedef std::deque<std::string> t_oneTraceback;
	typedef std::unordered_map<std::thread::id, t_oneTraceback> t_stckTraceback;
	t_stckTraceback stckTraceback;

	void PushTraceback(std::string&& strTrace);
	void PopTraceback();



	ParseInfo();
	~ParseInfo();

	NodeFunction* GetFunction(const t_string& strName);
};


struct ParseObj
{
	Lexer* pLexer;
	Node* pRoot;

	// only used during parsing/lexing for yyerror(), NOT during exec
	unsigned int iCurLine;
	t_string strCurFile;
};

#endif
