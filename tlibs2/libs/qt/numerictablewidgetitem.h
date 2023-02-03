/**
 * numeric table widget item
 * @author Tobias Weber <tweber@ill.fr>
 * @date Dec-2018
 * @license GPLv3, see 'LICENSE' file
 * @desc The present version was forked on 28-Dec-2018 and 21-Apr-2021 from my privately developed "misc" project (https://github.com/t-weber/misc).
 *
 * ----------------------------------------------------------------------------
 * tlibs
 * Copyright (C) 2017-2021  Tobias WEBER (Institut Laue-Langevin (ILL),
 *                          Grenoble, France).
 * Copyright (C) 2015-2017  Tobias WEBER (Technische Universitaet Muenchen
 *                          (TUM), Garching, Germany).
 * "misc" project
 * Copyright (C) 2017-2021  Tobias WEBER (privately developed).
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, version 3 of the License.
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

#ifndef __NUM_TABWIDGETITEM_H__
#define __NUM_TABWIDGETITEM_H__


#include <QtWidgets/QTableWidget>
#include "../str.h"


namespace tl2 {


template<class T = double>
class NumericTableWidgetItem : public QTableWidgetItem
{
public:
	/**
	 * construct a numeric table widget item from a numeric value
	 */
	explicit NumericTableWidgetItem(const T& val, std::streamsize prec = 6)
		: QTableWidgetItem(tl2::var_to_str(val, prec).c_str()),
		  m_prec{prec}
	{}


	/**
	 * construct a numeric table widget item from a string value
	 */
	explicit NumericTableWidgetItem(const QString& val, std::streamsize prec = 6)
		: QTableWidgetItem(val), m_prec{prec}
	{}


	/**
	 * construct a numeric table widget item from a string value
	 */
	explicit NumericTableWidgetItem(const std::string& val, std::streamsize prec = 6)
		: QTableWidgetItem(QString(val.c_str())), m_prec{prec}
	{}


	/**
	 * compare items
	 */
	virtual bool operator<(const QTableWidgetItem& item) const override
	{
		T val1 = GetValue();
		T val2 = tl2::str_to_var_parse<T>(item.text().toStdString());

		return val1 < val2;
	}


	/**
	 * get a cloned copy of this table item widget
	 */
	virtual QTableWidgetItem* clone() const override
	{
		auto item = new (std::nothrow)
			NumericTableWidgetItem<T>(GetValue(), m_prec);
		if(item)
			*static_cast<QTableWidgetItem*>(item) = *this;

		return item;
	};


	/**
	 * get the numerical value
	 */
	T GetValue() const
	{
		return tl2::str_to_var_parse<T>(text().toStdString());
	}


	/**
	 * set the numerical value
	 */
	void SetValue(T val)
	{
		QString str = tl2::var_to_str(val, m_prec).c_str();
		setText(str);
	}


private:
	// precision
	std::streamsize m_prec{6};
};


}
#endif
