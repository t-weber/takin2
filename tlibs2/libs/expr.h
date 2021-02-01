/**
 * tlibs2 -- simple LL(1) expression parser
 * @author Tobias Weber <tweber@ill.fr>
 * @date 28-mar-2020
 * @license GPLv3, see 'LICENSE' file
 * @desc The present version was forked on 28-Mar-2020 from my privately developed "misc" project (https://github.com/t-weber/misc).
 *
 * References:
 *   - https://de.wikipedia.org/wiki/LL(k)-Grammatik
 *   - https://www.cs.uaf.edu/~cs331/notes/FirstFollow.pdf
 */

#ifndef __TLIBS2_EXPR_PARSER_H__
#define __TLIBS2_EXPR_PARSER_H__

//#define TL2_USE_THREADS
#define TL2_USE_UNITS


#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <regex>
#include <unordered_map>
#include <vector>
#include <stack>
#include <memory>
#include <cmath>
#include <cstdint>

#ifdef TL2_USE_THREADS
	#include <future>
#endif

#if __has_include(<boost/math/special_functions/erf.hpp>)
	#include <boost/math/special_functions/erf.hpp>
	#define HAS_ERFINV
#endif

#ifdef TL2_USE_UNITS
	#include "units.h"
#endif


namespace tl2 {


template<class t_num=double>
t_num expr_modfunc(t_num t1, t_num t2)
{
	if constexpr(std::is_floating_point_v<t_num>)
		return std::fmod(t1, t2);
	else if constexpr(std::is_integral_v<t_num>)
		return t1%t2;
	else
		throw std::runtime_error{"Invalid type for mod function."};
}


template<class t_num=double>
t_num expr_binop(char op, t_num val_left, t_num val_right)
{
	switch(op)
	{
		case '+': return val_left + val_right;
		case '-': return val_left - val_right;
		case '*': return val_left * val_right;
		case '/': return val_left / val_right;
		case '%': return expr_modfunc<t_num>(val_left, val_right);
		case '^': return static_cast<t_num>(std::pow(val_left, val_right));
	}

	throw std::runtime_error{"Invalid binary operator."};
}


template<class t_num=double>
t_num expr_unop(char op, t_num val)
{
	switch(op)
	{
		case '+': return val;
		case '-': return -val;
	}

	throw std::runtime_error{"Invalid unary operator."};
}


template<typename t_num=double> class ExprParser;


// ------------------------------------------------------------------------
// vm
// ------------------------------------------------------------------------

template<typename t_num=double>
class ExprVM
{
public:
	enum class Op : std::uint8_t
	{
		NOP = 0,
		BINOP, UNOP,
		PUSH_VAR, PUSH_VAL,
		CALL,
	};


public:
	ExprVM(bool debug) : m_stack{}, m_debug{debug}
	{}


	t_num run(const std::uint8_t* code, std::size_t size, const ExprParser<t_num>& context)
	{
		for(std::size_t ip=0; ip<size;)
		{
			const Op op = *reinterpret_cast<const Op*>(code + ip++);

			switch(op)
			{
				case Op::NOP:
				{
					break;
				}
				case Op::BINOP:
				{
					char binop = *reinterpret_cast<const char*>(code + ip++);

					t_num val2 = m_stack.top(); m_stack.pop();
					t_num val1 = m_stack.top(); m_stack.pop();

					t_num result = expr_binop<t_num>(binop, val1, val2);
					m_stack.push(result);

					if(m_debug)
						std::cout << val1 << " " << binop << " " << val2 << " = " << result << std::endl;
					break;
				}
				case Op::UNOP:
				{
					char unop = *reinterpret_cast<const char*>(code + ip++);

					t_num val = m_stack.top(); m_stack.pop();
					t_num result = expr_unop<t_num>(unop, val);
					m_stack.push(result);

					if(m_debug)
						std::cout << unop << " " << val << " = " << result << std::endl;

					break;
				}
				case Op::PUSH_VAR:
				{
					std::size_t len = *reinterpret_cast<const std::size_t*>(code + ip);
					ip += sizeof(len);

					std::string var(reinterpret_cast<const char*>(code + ip),
						reinterpret_cast<const char*>(code + ip+len));
					ip += len;

					t_num val = context.get_var_or_const(var);
					m_stack.push(val);

					if(m_debug)
						std::cout << var << " = " << val << std::endl;

					break;
				}
				case Op::PUSH_VAL:
				{
					t_num val = *reinterpret_cast<const t_num*>(code + ip);
					ip += sizeof(val);

					m_stack.push(val);


					if(m_debug)
						std::cout << val << " = " << val << std::endl;
					break;
				}
				case Op::CALL:
				{
					std::uint8_t numargs = *reinterpret_cast<const std::uint8_t*>(code + ip);
					ip += sizeof(numargs);

					std::size_t len = *reinterpret_cast<const std::size_t*>(code + ip);
					ip += sizeof(len);

					std::string fkt(reinterpret_cast<const char*>(code + ip),
						reinterpret_cast<const char*>(code + ip+len));
					ip += len;

					if(numargs == 0)
					{
						t_num result = context.call_func0(fkt);
						m_stack.push(result);

						if(m_debug)
							std::cout << fkt << "() = " << result << std::endl;
					}
					else if(numargs == 1)
					{
						t_num argval = m_stack.top(); m_stack.pop();
						t_num result = context.call_func1(fkt, argval);

						m_stack.push(result);

						if(m_debug)
							std::cout << fkt << "(" << argval << ") = " << result << std::endl;
					}
					else if(numargs == 2)
					{
						t_num arg2val = m_stack.top(); m_stack.pop();
						t_num arg1val = m_stack.top(); m_stack.pop();

						t_num result = context.call_func2(fkt, arg1val, arg2val);

						m_stack.push(result);

						if(m_debug)
							std::cout << fkt << "(" << arg1val << ", " << arg2val << ") = " << result << std::endl;
					}
					else
					{
						throw std::runtime_error("Invalid function call.");
					}

					break;
				}
				default:
				{
					throw std::runtime_error("Invalid opcode.");
				}
			}
		}

		if(m_stack.size() != 1)
			throw std::runtime_error("Result not on stack");
		t_num result = m_stack.top(); m_stack.pop();
		return result;
	}


private:
	std::stack<t_num> m_stack;
	bool m_debug = false;
};

// ------------------------------------------------------------------------



// ------------------------------------------------------------------------
// ast
// ------------------------------------------------------------------------
template<typename t_num=double>
class ExprAST
{
public:
	virtual ~ExprAST() = default;

	virtual t_num eval(const ExprParser<t_num>&) const = 0;
	virtual void codegen(std::ostream& ostr) const = 0;

	virtual void print(std::ostream& ostr=std::cout, std::size_t indent=0) const = 0;
};


template<typename t_num=double>
class ExprASTBinOp : public ExprAST<t_num>
{
public:
	ExprASTBinOp(char op, const std::shared_ptr<ExprAST<t_num>>& left, const std::shared_ptr<ExprAST<t_num>>& right)
		: m_op{op}, m_left{left}, m_right{right}
	{}

	ExprASTBinOp(const ExprASTBinOp<t_num>&) = delete;
	const ExprASTBinOp<t_num>& operator=(const ExprASTBinOp<t_num>&) = delete;

	virtual ~ExprASTBinOp() = default;

	virtual t_num eval(const ExprParser<t_num>& context) const override
	{
#ifdef TL2_USE_THREADS
		auto fut_left = std::async([this, &context]() { return m_left->eval(context); });
		const t_num val_right = m_right->eval(context);
		const t_num val_left = fut_left.get();
#else
		const t_num val_left = m_left->eval(context);
		const t_num val_right = m_right->eval(context);
#endif

		return expr_binop<t_num>(m_op, val_left, val_right);
	}

	virtual void codegen(std::ostream& ostr) const override
	{
		m_left->codegen(ostr);
		m_right->codegen(ostr);

		auto op = ExprVM<t_num>::Op::BINOP;
		ostr.write(reinterpret_cast<const char*>(&op), sizeof(op));
		ostr.write(&m_op, sizeof(m_op));
	}

	virtual void print(std::ostream& ostr=std::cout, std::size_t indent=0) const override
	{
		for(std::size_t i=0; i<indent; ++i) ostr << " | ";
		ostr << "binary operator " << m_op << "\n";
		m_left->print(ostr, indent+1);
		m_right->print(ostr, indent+1);
	}

private:
	char m_op{'+'};
	std::shared_ptr<ExprAST<t_num>> m_left{}, m_right{};
};


template<typename t_num=double>
class ExprASTUnOp : public ExprAST<t_num>
{
public:
	ExprASTUnOp(char op, const std::shared_ptr<ExprAST<t_num>>& child)
		: m_op{op}, m_child{child}
	{}

	ExprASTUnOp(const ExprASTUnOp<t_num>&) = delete;
	const ExprASTUnOp<t_num>& operator=(const ExprASTUnOp<t_num>&) = delete;

	virtual ~ExprASTUnOp() = default;

	virtual t_num eval(const ExprParser<t_num>& context) const override
	{
		const t_num val = m_child->eval(context);
		return expr_unop<t_num>(m_op, val);
	}

	virtual void codegen(std::ostream& ostr) const override
	{
		m_child->codegen(ostr);

		auto op = ExprVM<t_num>::Op::UNOP;
		ostr.write(reinterpret_cast<const char*>(&op), sizeof(op));
		ostr.write(&m_op, sizeof(m_op));
	}

	virtual void print(std::ostream& ostr=std::cout, std::size_t indent=0) const override
	{
		for(std::size_t i=0; i<indent; ++i) ostr << " | ";
		ostr << "unary operator " << m_op << "\n";
		m_child->print(ostr, indent+1);
	}

private:
	char m_op{'+'};
	std::shared_ptr<ExprAST<t_num>> m_child{};
};


template<typename t_num=double>
class ExprASTVar : public ExprAST<t_num>
{
public:
	ExprASTVar(const std::string& name) : m_name{name}
	{}

	ExprASTVar(const ExprASTVar<t_num>&) = delete;
	const ExprASTVar<t_num>& operator=(const ExprASTVar<t_num>&) = delete;

	virtual ~ExprASTVar() = default;

	virtual t_num eval(const ExprParser<t_num>& context) const override
	{
		return context.get_var_or_const(m_name);
	}

	virtual void codegen(std::ostream& ostr) const override
	{
		auto op = ExprVM<t_num>::Op::PUSH_VAR;
		std::size_t len = m_name.length();

		ostr.write(reinterpret_cast<const char*>(&op), sizeof(op));
		ostr.write(reinterpret_cast<const char*>(&len), sizeof(len));
		ostr.write(m_name.c_str(), m_name.length());
	}

	virtual void print(std::ostream& ostr=std::cout, std::size_t indent=0) const override
	{
		for(std::size_t i=0; i<indent; ++i) ostr << " | ";
		ostr << "variable \"" << m_name << "\"\n";
	}

private:
	std::string m_name{};
};


template<typename t_num=double>
class ExprASTValue : public ExprAST<t_num>
{
public:
	ExprASTValue(t_num val) : m_val{val}
	{}

	ExprASTValue(const ExprASTValue<t_num>&) = delete;
	const ExprASTValue<t_num>& operator=(const ExprASTValue<t_num>&) = delete;

	virtual ~ExprASTValue() = default;

	virtual t_num eval(const ExprParser<t_num>&) const override
	{
		return m_val;
	}

	virtual void codegen(std::ostream& ostr) const override
	{
		auto op = ExprVM<t_num>::Op::PUSH_VAL;

		ostr.write(reinterpret_cast<const char*>(&op), sizeof(op));
		ostr.write(reinterpret_cast<const char*>(&m_val), sizeof(m_val));
	}

	virtual void print(std::ostream& ostr=std::cout, std::size_t indent=0) const override
	{
		for(std::size_t i=0; i<indent; ++i) ostr << " | ";
		ostr << "value " << m_val << "\n";
	}

private:
	t_num m_val{};
};


template<typename t_num=double>
class ExprASTCall : public ExprAST<t_num>
{
public:
	ExprASTCall(const std::string& name)
		: m_name{name}
	{}

	ExprASTCall(const std::string& name, const std::shared_ptr<ExprAST<t_num>>& arg1)
		: m_name{name}
	{
		m_args.push_back(arg1);
	}

	ExprASTCall(const std::string& name,
		const std::shared_ptr<ExprAST<t_num>>& arg1, const std::shared_ptr<ExprAST<t_num>>& arg2)
		: m_name{name}
	{
		m_args.push_back(arg1);
		m_args.push_back(arg2);
	}

	ExprASTCall(const ExprASTCall<t_num>&) = delete;
	const ExprASTCall<t_num>& operator=(const ExprASTCall<t_num>&) = delete;

	virtual ~ExprASTCall() = default;

	virtual t_num eval(const ExprParser<t_num>& context) const override
	{
		if(m_args.size() == 0)
		{
			return context.call_func0(m_name);
		}
		else if(m_args.size() == 1)
		{
			return context.call_func1(m_name, m_args[0]->eval(context));
		}
		else if(m_args.size() == 2)
		{
#ifdef TL2_USE_THREADS
			auto fut_arg0 = std::async([this, &context]() { return m_args[0]->eval(context); });
			const t_num arg1 = m_args[1]->eval(context);
			const t_num arg0 = fut_arg0.get();
#else
			const t_num arg0 = m_args[0]->eval(context);
			const t_num arg1 = m_args[1]->eval(context);
#endif

			return context.call_func2(m_name, arg0, arg1);
		}

		throw std::runtime_error("Invalid function call.");
	}

	virtual void codegen(std::ostream& ostr) const override
	{
		for(const auto& arg : m_args)
			arg->codegen(ostr);

		auto op = ExprVM<t_num>::Op::CALL;
		std::uint8_t numargs = static_cast<std::uint8_t>(m_args.size());
		std::size_t len = m_name.length();

		ostr.write(reinterpret_cast<const char*>(&op), sizeof(op));
		ostr.write(reinterpret_cast<const char*>(&numargs), sizeof(numargs));
		ostr.write(reinterpret_cast<const char*>(&len), sizeof(len));
		ostr.write(m_name.c_str(), len);
	}

	virtual void print(std::ostream& ostr=std::cout, std::size_t indent=0) const override
	{
		for(std::size_t i=0; i<indent; ++i) ostr << " | ";
		ostr << "function call \"" << m_name << "\"\n";

		for(const auto& arg : m_args)
			arg->print(ostr, indent+1);
	}

private:
	std::string m_name{};
	std::vector<std::shared_ptr<ExprAST<t_num>>> m_args{};
};


// ------------------------------------------------------------------------



template<typename t_num>
class ExprParser
{
	friend class ExprASTCall<t_num>;
	friend class ExprASTVar<t_num>;
	friend class ExprVM<t_num>;


public:
	ExprParser(bool debug=false)
		: m_debug{debug}, m_ast{}, m_code{},
		m_vars{}, m_consts{}, m_funcs0{}, m_funcs1{}, m_funcs2{},
		m_istr{}, m_lookahead_text{}
	{
		register_funcs();
		register_consts();
	}


	bool parse(const std::string& str, bool codegen=true)
	{
		m_code.clear();

		m_istr = std::make_shared<std::istringstream>(str);
		next_lookahead();
		m_ast = plus_term();

		// check if there would be are more tokens available?
		next_lookahead();
		bool at_eof = (m_lookahead == (int)Token::TOK_INVALID || m_lookahead == (int)Token::TOK_END);
		if(!at_eof)
			throw std::underflow_error("Not all input tokens have been consumed.");

		bool ok = !!m_ast;
		if(ok)
		{
			if(m_debug)
			{
				m_ast->print(std::cout);
				std::cout << std::endl;
			}

			if(codegen)
			{
				std::stringstream code;
				m_ast->codegen(code);

				code.seekg(0, std::stringstream::end);
				auto code_size = code.tellg();
				code.seekg(0, std::stringstream::beg);

				m_code.resize(code_size);
				code.read(reinterpret_cast<char*>(m_code.data()), code_size);
			}
		}

		return ok;
	}


	t_num eval()
	{
		// is compiled code available?
		if(m_code.size())
		{
			// run generated code
			ExprVM<t_num> vm{m_debug};
			return vm.run(m_code.data(), m_code.size(), *this);
		}

		if(m_debug)
			std::cerr << "Warning: No code available, interpreting AST." << std::endl;
		if(!m_ast)
			throw std::runtime_error("Invalid AST.");
		return m_ast->eval(*this);
	}


protected:
	// ------------------------------------------------------------------------
	// tables / functions
	// ------------------------------------------------------------------------

	void register_funcs()
	{
		// common functions
		register_func1("abs", std::abs);
		register_func2("mod", expr_modfunc<t_num>);

		// real functions
		if constexpr(std::is_floating_point_v<t_num>)
		{
			register_func1("sin", std::sin);
			register_func1("cos", std::cos);
			register_func1("tan", std::tan);
			register_func1("asin", std::asin);
			register_func1("acos", std::acos);
			register_func1("atan", std::atan);
			register_func1("sinh", std::sinh);
			register_func1("cosh", std::cosh);
			register_func1("tanh", std::tanh);
			register_func1("asinh", std::asinh);
			register_func1("acosh", std::acosh);
			register_func1("atanh", std::atanh);
			register_func1("sqrt", std::sqrt);
			register_func1("cbrt", std::cbrt);
			register_func1("exp", std::exp);
			register_func1("log", std::log);
			register_func1("log2", std::log2);
			register_func1("log10", std::log10);
			register_func1("erf", std::erf);
			register_func1("erfc", std::erfc);
			register_func1("round", std::round);
			register_func1("ceil", std::ceil);
			register_func1("floor", std::floor);
#ifdef HAS_ERFINV
			register_func1("erf_inv", boost::math::erf_inv);
#endif

			register_func2("pow", std::pow);
			register_func2("atan2", std::atan2);
		}

		// integer functions
		else if constexpr(std::is_integral_v<t_num>)
		{
			register_func2("pow", [](t_num t1, t_num t2) -> t_num { return t_num(std::pow(t1, t2)); } );
		}
	}


	// call function with zero parameters
	t_num call_func0(const std::string& strName) const
	{
		return m_funcs0.at(strName)();
	}


	// call function with one parameter
	t_num call_func1(const std::string& strName, t_num t) const
	{
		return m_funcs1.at(strName)(t);
	}


	// call function with two parameters
	t_num call_func2(const std::string& strName, t_num t1, t_num t2) const
	{
		return m_funcs2.at(strName)(t1, t2);
	}


	// register constants
	void register_consts()
	{
		// real constants
		if constexpr(std::is_floating_point_v<t_num>)
		{
#ifdef TL2_USE_UNITS
			register_const("pi", __pi<t_num>);
			register_const("hbar",  t_num(hbar<t_num>/meV<t_num>/sec<t_num>));	// hbar in [meV s]
			register_const("kB",  t_num(kB<t_num>/meV<t_num>*kelvin<t_num>));		// kB in [meV / K]
#endif
		}

		// integer constants
		else if constexpr(std::is_integral_v<t_num>)
		{
		}
	}


	// get variable or constant
	t_num get_var_or_const(const std::string& strName) const
	{
		// look for variable
		if(const auto iterVar = m_vars.find(strName); iterVar != m_vars.end())
			return iterVar->second;

		// else look for constant
		if(const auto iterConst = m_consts.find(strName); iterConst != m_consts.end())
			return iterConst->second;

		throw std::runtime_error("Variable or constant \"" + strName + "\" not found.");
	}


	// ------------------------------------------------------------------------



	// ------------------------------------------------------------------------
	// Lexer
	// ------------------------------------------------------------------------
	enum class Token : int
	{
		TOK_NUM		= 1000,
		TOK_IDENT	= 1001,
		TOK_END		= 1002,

		TOK_INVALID	= 10000,
	};


	/**
	 * find all matching tokens for input string
	 */
	std::vector<std::pair<int, t_num>> get_matching_tokens(const std::string& str)
	{
		std::vector<std::pair<int, t_num>> matches;

		if constexpr(std::is_floating_point_v<t_num>)
		{	// real
			std::regex regex{"[0-9]+(\\.[0-9]*)?|\\.[0-9]+([eE][-+]?[0-9]+)?"};
			std::smatch smatch;
			if(std::regex_match(str, smatch, regex))
			{
				t_num val{};
				std::istringstream{str} >> val;
				matches.push_back(std::make_pair((int)Token::TOK_NUM, val));
			}
		}
		else if constexpr(std::is_integral_v<t_num>)
		{	// real
			std::regex regex{"[0-9]+"};
			std::smatch smatch;
			if(std::regex_match(str, smatch, regex))
			{
				t_num val{};
				std::istringstream{str} >> val;
				matches.push_back(std::make_pair((int)Token::TOK_NUM, val));
			}
		}
		else
		{
			throw std::invalid_argument("Invalid number type.");
		}

		{	// ident
			std::regex regex{"[A-Za-z_][A-Za-z0-9_]*"};
			std::smatch smatch;
			if(std::regex_match(str, smatch, regex))
				matches.push_back(std::make_pair((int)Token::TOK_IDENT, 0.));
		}

		{	// tokens represented by themselves
			if(str == "+" || str == "-" || str == "*" || str == "/" ||
				str == "%" || str == "^" || str == "(" || str == ")" || str == ",")
				matches.push_back(std::make_pair((int)str[0], 0.));
		}

		return matches;
	}


	/**
	 * @return [token, yylval, yytext]
	 */
	std::tuple<int, t_num, std::string> lex()
	{
		std::string input, longest_input;
		std::vector<std::pair<int, t_num>> longest_matching;

		// find longest matching token
		while(1)
		{
			char c = m_istr->get();

			if(m_istr->eof())
				break;
			// if outside any other match...
			if(longest_matching.size() == 0)
			{
				// ...ignore white spaces
				if(c==' ' || c=='\t')
					continue;
				// ...end on new line
				if(c=='\n')
					return std::make_tuple((int)Token::TOK_END, t_num{0}, longest_input);
			}

			input += c;
			auto matching = get_matching_tokens(input);
			if(matching.size())
			{
				longest_input = input;
				longest_matching = matching;

				if(m_istr->peek() == std::char_traits<char>::eof())
					break;
			}
			else
			{
				// no more matches
				m_istr->putback(c);
				break;
			}
		}

		// at EOF
		if(longest_matching.size() == 0 && (input.length() == 0 || m_istr->eof()))
		{
			return std::make_tuple((int)Token::TOK_END, t_num{0}, longest_input);
		}

		// nothing matches
		if(longest_matching.size() == 0)
		{
			std::ostringstream ostr;
			ostr << "Invalid input in lexer: \"" << input << "\".";
			throw std::runtime_error(ostr.str());
		}

		// several possible matches
		if(longest_matching.size() > 1)
		{
			std::ostringstream ostr;
			ostr << "Warning: Ambiguous match in lexer for token \"" << longest_input << "\".";
			throw std::runtime_error(ostr.str());
		}

		// found match
		return std::make_tuple((int)std::get<0>(longest_matching[0]), std::get<1>(longest_matching[0]), longest_input);
	}
	// ------------------------------------------------------------------------



	// ----------------------------------------------------------------------------
	// Lexer interface
	// ----------------------------------------------------------------------------
	void next_lookahead()
	{
		std::tie(m_lookahead, m_lookahead_val, m_lookahead_text) = lex();
	}


	void match(int expected)
	{
		if(m_lookahead != expected)
		{
			std::ostringstream ostr;
			ostr << "Could not match symbol! Expected: " << expected << ", got: " << m_lookahead << ".";
			throw std::runtime_error(ostr.str());
		}
	}
	// ----------------------------------------------------------------------------



	// ----------------------------------------------------------------------------
	// Productions
	// ----------------------------------------------------------------------------
	/**
	 * +,- terms
	 * (lowest precedence, 1)
	 */
	std::shared_ptr<ExprAST<t_num>> plus_term()
	{
		// plus_term -> mul_term plus_term_rest
		if(m_lookahead == '(' || m_lookahead == (int)Token::TOK_NUM || m_lookahead == (int)Token::TOK_IDENT)
		{
			auto term_val = mul_term();
			auto expr_rest_val = plus_term_rest(term_val);

			return expr_rest_val;
		}
		else if(m_lookahead == '+')	// unary +
		{
			next_lookahead();
			auto term_val = mul_term();
			auto expr_rest_val = plus_term_rest(term_val);

			return expr_rest_val;
		}
		else if(m_lookahead == '-')	// unary -
		{
			next_lookahead();
			auto term_val = std::make_shared<ExprASTUnOp<t_num>>('-', mul_term());
			auto expr_rest_val = plus_term_rest(term_val);

			return expr_rest_val;
		}

		std::ostringstream ostr;
		if(m_lookahead == 0 || m_lookahead == EOF)
			ostr << "EOF in " << __func__ << ".";
		else
			ostr << "Invalid lookahead in " << __func__ << ": " << m_lookahead << ".";
		throw std::runtime_error(ostr.str());
	}


	std::shared_ptr<ExprAST<t_num>> plus_term_rest(const std::shared_ptr<ExprAST<t_num>>& arg)
	{
		// plus_term_rest -> '+' mul_term plus_term_rest
		if(m_lookahead == '+')
		{
			next_lookahead();
			auto term_val = std::make_shared<ExprASTBinOp<t_num>>('+', arg, mul_term());
			auto expr_rest_val = plus_term_rest(term_val);

			return expr_rest_val;
		}

		// plus_term_rest -> '-' mul_term plus_term_rest
		else if(m_lookahead == '-')
		{
			next_lookahead();
			auto term_val = std::make_shared<ExprASTBinOp<t_num>>('-', arg, mul_term());
			auto expr_rest_val = plus_term_rest(term_val);

			return expr_rest_val;
		}
		// plus_term_rest -> epsilon
		else if(m_lookahead == ')' || m_lookahead == (int)Token::TOK_END || m_lookahead == ',')
		{
			return arg;
		}

		std::ostringstream ostr;
		ostr << "Invalid lookahead in " << __func__ << ": " << m_lookahead << ".";
		throw std::runtime_error(ostr.str());
	}


	/**
	 * *,/,% terms
	 * (precedence 2)
	 */
	std::shared_ptr<ExprAST<t_num>> mul_term()
	{
		// mul_term -> pow_term mul_term_rest
		if(m_lookahead == '(' || m_lookahead == (int)Token::TOK_NUM || m_lookahead == (int)Token::TOK_IDENT)
		{
			auto factor_val = pow_term();
			auto term_rest_val = mul_term_rest(factor_val);

			return term_rest_val;
		}

		std::ostringstream ostr;
		ostr << "Invalid lookahead in " << __func__ << ": " << m_lookahead << ".";
		throw std::runtime_error(ostr.str());
	}


	std::shared_ptr<ExprAST<t_num>> mul_term_rest(const std::shared_ptr<ExprAST<t_num>>& arg)
	{
		// mul_term_rest -> '*' pow_term mul_term_rest
		if(m_lookahead == '*')
		{
			next_lookahead();
			auto factor_val = std::make_shared<ExprASTBinOp<t_num>>('*', arg, pow_term());
			auto term_rest_val = mul_term_rest(factor_val);

			return term_rest_val;
		}

		// mul_term_rest -> '/' pow_term mul_term_rest
		else if(m_lookahead == '/')
		{
			next_lookahead();
			auto factor_val = std::make_shared<ExprASTBinOp<t_num>>('/', arg, pow_term());
			auto term_rest_val = mul_term_rest(factor_val);

			return term_rest_val;
		}

		// mul_term_rest -> '%' pow_term mul_term_rest
		else if(m_lookahead == '%')
		{
			next_lookahead();
			auto factor_val = std::make_shared<ExprASTBinOp<t_num>>('%', arg, pow_term());
			auto term_rest_val = mul_term_rest(factor_val);

			return term_rest_val;
		}

		// mul_term_rest -> epsilon
		else if(m_lookahead == '+' || m_lookahead == '-' || m_lookahead == ')'
			|| m_lookahead == (int)Token::TOK_END || m_lookahead == ',')
		{
			return arg;
		}

		std::ostringstream ostr;
		ostr << "Invalid lookahead in " << __func__ << ": " << m_lookahead << ".";
		throw std::runtime_error(ostr.str());
	}


	/**
	 * ^ terms
	 * (precedence 3)
	 */
	std::shared_ptr<ExprAST<t_num>> pow_term()
	{
		// pow_term -> factor pow_term_rest
		if(m_lookahead == '(' || m_lookahead == (int)Token::TOK_NUM || m_lookahead == (int)Token::TOK_IDENT)
		{
			auto factor_val = factor();
			auto term_rest_val = pow_term_rest(factor_val);

			return term_rest_val;
		}

		std::ostringstream ostr;
		ostr << "Invalid lookahead in " << __func__ << ": " << m_lookahead << ".";
		throw std::runtime_error(ostr.str());
	}


	std::shared_ptr<ExprAST<t_num>> pow_term_rest(const std::shared_ptr<ExprAST<t_num>>& arg)
	{
		// pow_term_rest -> '^' factor pow_term_rest
		if(m_lookahead == '^')
		{
			next_lookahead();
			auto factor_val = std::make_shared<ExprASTBinOp<t_num>>('^', arg, factor());
			auto term_rest_val = pow_term_rest(factor_val);

			return term_rest_val;
		}

		// pow_term_rest -> epsilon
		else if(m_lookahead == '+' || m_lookahead == '-' || m_lookahead == ')'
			|| m_lookahead == (int)Token::TOK_END || m_lookahead == ','
			|| m_lookahead == '*' || m_lookahead == '/' || m_lookahead == '%')
		{
			return arg;
		}

		std::ostringstream ostr;
		ostr << "Invalid lookahead in " << __func__ << ": " << m_lookahead << ".";
		throw std::runtime_error(ostr.str());
	}


	/**
	 * () terms, real factor or identifier
	 * (highest precedence, 4)
	 */
	std::shared_ptr<ExprAST<t_num>> factor()
	{
		// factor -> '(' plus_term ')'
		if(m_lookahead == '(')
		{
			next_lookahead();
			auto expr_val = plus_term();
			match(')');
			next_lookahead();

			return expr_val;
		}

		// factor -> TOK_NUM
		else if(m_lookahead == (int)Token::TOK_NUM)
		{
			t_num val = m_lookahead_val;
			next_lookahead();

			return std::make_shared<ExprASTValue<t_num>>(val);
		}

		// factor -> TOK_IDENT
		else if(m_lookahead == (int)Token::TOK_IDENT)
		{
			const std::string ident = m_lookahead_text;
			next_lookahead();

			// function call
			// using next m_lookahead, grammar still ll(1)?
			if(m_lookahead == '(')
			{
				next_lookahead();

				// 0-argument function
				// factor -> TOK_IDENT '(' ')'
				if(m_lookahead == ')')
				{
					next_lookahead();

					return std::make_shared<ExprASTCall<t_num>>(ident);
				}

				// function with arguments
				else
				{
					// first argument
					auto expr_val1 = plus_term();

					// one-argument-function
					// factor -> TOK_IDENT '(' plus_term ')'
					if(m_lookahead == ')')
					{
						next_lookahead();

						return std::make_shared<ExprASTCall<t_num>>(ident, expr_val1);
					}

					// two-argument-function
					// factor -> TOK_IDENT '(' plus_term ',' plus_term ')'
					else if(m_lookahead == ',')
					{
						next_lookahead();
						auto expr_val2 = plus_term();
						match(')');
						next_lookahead();

						return std::make_shared<ExprASTCall<t_num>>(ident, expr_val1, expr_val2);
					}
					else
					{
						std::ostringstream ostr;
						ostr << "Invalid function call to \"" << ident << "\".";
						throw std::runtime_error(ostr.str());
					}
				}
			}

			// variable lookup
			else
			{
				// register the variable if it doesn't yet exist
				if(m_vars.find(ident) == m_vars.end() && m_consts.find(ident) == m_consts.end())
					register_var(ident, t_num{});
				return std::make_shared<ExprASTVar<t_num>>(ident);
			}
		}

		std::ostringstream ostr;
		ostr << "Invalid lookahead in " << __func__ << ": " << m_lookahead << ".";
		throw std::runtime_error(ostr.str());
	}
	// ----------------------------------------------------------------------------


public:
	// register a function with no parameters
	void register_func0(const std::string& name, t_num(*fkt)())
	{
		m_funcs0.emplace(std::make_pair(name, static_cast<t_num(*)()>(fkt)));
	}


	// register a function with one parameter
	void register_func1(const std::string& name, t_num(*fkt)(t_num))
	{
		m_funcs1.emplace(std::make_pair(name, static_cast<t_num(*)(t_num)>(fkt)));
	}


	// register a function with two parameters
	void register_func2(const std::string& name, t_num(*fkt)(t_num, t_num))
	{
		m_funcs2.emplace(std::make_pair(name, static_cast<t_num(*)(t_num, t_num)>(fkt)));
	}


	// register a variable
	void register_var(const std::string& name, t_num val)
	{
		// overwrite value if key already exists
		if(auto [iter, ok] = m_vars.emplace(std::make_pair(name, val)); !ok)
			iter->second = val;
	}


	// register a constant
	void register_const(const std::string& name, t_num val)
	{
		// overwrite value if key already exists
		if(auto [iter, ok] = m_consts.emplace(std::make_pair(name, val)); !ok)
			iter->second = val;
	}


	const std::unordered_map<std::string, t_num>& get_vars() const
	{
		return m_vars;
	}


private:
	bool m_debug = false;

	// ast root
	std::shared_ptr<ExprAST<t_num>> m_ast{};

	// generated code
	std::vector<std::uint8_t> m_code{};

	// variables and constants
	std::unordered_map<std::string, t_num> m_vars{}, m_consts{};

	// functions
	std::unordered_map<std::string, t_num(*)()> m_funcs0{};
	std::unordered_map<std::string, t_num(*)(t_num)> m_funcs1{};
	std::unordered_map<std::string, t_num(*)(t_num, t_num)> m_funcs2{};


private:
	std::shared_ptr<std::istream> m_istr{};

	int m_lookahead = (int)Token::TOK_INVALID;
	t_num m_lookahead_val = 0;
	std::string m_lookahead_text = "";
};


}
#endif
