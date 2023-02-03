/**
 * Cache dialog
 * @author Tobias Weber <tobias.weber@tum.de>
 * @date 21-oct-2014
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

#ifndef __NET_CACHE_DLG_H__
#define __NET_CACHE_DLG_H__

#include <map>
#include <QDialog>
#include <QSettings>
#include <QTimer>
#include "ui/ui_netcache.h"
#include "libs/globals.h"


enum class CacheValType : int
{
	UNKNOWN,

	TIMER,
	PRESET,
	COUNTER,

	LIVE_PLOT,
};

struct CacheVal
{
	std::string strVal;
	t_real_glob dTimestamp = t_real_glob(-1);

	CacheValType ty = CacheValType::UNKNOWN;
};

typedef std::map<std::string, CacheVal> t_mapCacheVal;


class NetCacheDlg : public QDialog, Ui::NetCacheDlg
{ Q_OBJECT
protected:
	QSettings *m_pSettings = 0;

	static const int s_iTimer = 1000;
	QTimer m_timer;

protected:
	virtual void hideEvent(QHideEvent *pEvt) override;
	virtual void showEvent(QShowEvent *pEvt) override;
	virtual void accept() override;

	void UpdateAge(int iRow=-1);

public:
	NetCacheDlg(QWidget* pParent=0, QSettings* pSett=0);
	virtual ~NetCacheDlg();

protected slots:
	void UpdateTimer();

public slots:
	void ClearAll();
	void UpdateValue(const std::string& strKey, const CacheVal& val);
	void UpdateAll(const t_mapCacheVal& map);

//signals:
//	void UpdatedValue(const std::string& strKey, const CacheVal& val);
};

#endif
