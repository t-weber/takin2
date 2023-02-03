/**
 * Interface for connections to instruments
 * @author Tobias Weber <tobias.weber@tum.de>
 * @date mar-2016
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
 
#ifndef __TAKIN_NET_IF_H__
#define __TAKIN_NET_IF_H__

#include <string>
#include <QObject>
#include <QString>

#include "tasoptions.h"
#include "dialogs/NetCacheDlg.h"


class NetCache : public QObject
{ Q_OBJECT
	public:
		virtual ~NetCache() {};

		virtual void connect(const std::string& strHost, const std::string& strPort, 
			const std::string& strUser, const std::string& strPass) = 0;
		virtual void disconnect() = 0;

		virtual void refresh() = 0;

	signals:
		void connected(const QString& strHost, const QString& strSrv);
		void disconnected();

		void vars_changed(const CrystalOptions& crys, const TriangleOptions& triag);

		void updated_cache_value(const std::string& strKey, const CacheVal& val);
		void cleared_cache();
};

#endif
