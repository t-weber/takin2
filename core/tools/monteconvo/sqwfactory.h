/**
 * factory and plugin interface for S(q,w) models
 * @author Tobias Weber <tobias.weber@tum.de>
 * @date 2016 -- 2018
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

#ifndef __SQW_FACT_H__
#define __SQW_FACT_H__

#include "sqwbase.h"
#include "tlibs/file/prop.h"

#include <tuple>
#include <memory>
#include <vector>
#include <string>


extern std::shared_ptr<SqwBase> construct_sqw(const std::string& strName,
	const std::string& strConfigFile);

// [identifier, long name]
extern std::vector<std::tuple<std::string, std::string>> get_sqw_names();


extern void unload_sqw_plugins();
extern void load_sqw_plugins();


// ----------------------------------------------------------------------------
// saving and loading of parameters

extern bool save_sqw_params(const SqwBase* pSqw,
	std::map<std::string, std::string>& mapConf, const std::string& strXmlRoot);

extern bool load_sqw_params(SqwBase* pSqw,
	const tl::Prop<std::string>& xml, const std::string& strXmlRoot);
// ----------------------------------------------------------------------------

#endif
