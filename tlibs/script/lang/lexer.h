/**
 * Lexer
 * @author Tobias Weber <tobias.weber@tum.de>
 * @date 2013
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

#ifndef __MIEZE_LEXER__
#define __MIEZE_LEXER__

#include "types.h"

#include <istream>
#include <string>
#include <vector>

enum TokenType : unsigned int
{
	LEX_TOKEN_INVALID,
	LEX_TOKEN_END,

	LEX_TOKEN_DOUBLE,
	LEX_TOKEN_INT,
	LEX_TOKEN_STRING,
	LEX_TOKEN_IDENT,
	LEX_TOKEN_CHAROP,

	LEX_TOKEN_IF,
	LEX_TOKEN_ELSE,
	LEX_TOKEN_FOR,
	LEX_TOKEN_WHILE,

	LEX_TOKEN_RETURN,
	LEX_TOKEN_BREAK,
	LEX_TOKEN_CONTINUE,

	LEX_TOKEN_LOG_AND,
	LEX_TOKEN_LOG_OR,
	LEX_TOKEN_LOG_NOT,
	LEX_TOKEN_LOG_EQ,
	LEX_TOKEN_LOG_NEQ,
	LEX_TOKEN_LOG_LESS,
	LEX_TOKEN_LOG_GREATER,
	LEX_TOKEN_LOG_LEQ,
	LEX_TOKEN_LOG_GEQ,

	LEX_TOKEN_GLOBAL
};

struct Token
{
	TokenType type;

	t_char cOp;
	t_int iVal;
	t_real dVal;
	t_string strVal;

	unsigned int iLine;

	Token() : type(LEX_TOKEN_INVALID), cOp(0), iVal(0), dVal(0), iLine(0)
	{}
};

class Lexer
{
protected:
	bool m_bOk;
	t_string m_strWhitespace, m_strSep;

	unsigned int m_iLexPos;
	unsigned int m_iNumToks;
	std::vector<Token> m_vecToks;
	Token m_tokEnd;

	t_string m_strFile;

	void FixTokens();

public:
	Lexer();
	Lexer(const t_string& str, const t_char* pcFile=0);
	virtual ~Lexer();

	void load(const t_string& strInput);
	void print();
	const Token& lex();

	unsigned int GetNumTokens() const { return m_vecToks.size(); }
	const Token& GetToken(unsigned int i) const { return m_vecToks[i]; }

	static t_string RemoveComments(const t_string& strInput);
	static std::vector<t_string> GetStringTable(const t_string& strInput);
	static void ReplaceEscapes(t_string& str);

	bool IsOk() const { return m_bOk; }
};

#endif
