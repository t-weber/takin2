/**
 * TcpClient
 * @author Tobias Weber <tobias.weber@tum.de>
 * @date aug-2014
 * @license GPLv2 or GPLv3
 *
 * @desc based on Boost's example chat client (c) 2003-2015 by C. M. Kohlhoff:
 * @desc http://www.boost.org/doc/libs/1_61_0/doc/html/boost_asio/example/cpp11/chat/chat_client.cpp
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

#include "tcp.h"
#include "tcp_impl.h"

template class tl::TcpTxtClient<>;
template class tl::TcpTxtServer<>;
