/**
 * Command line
 * @author Tobias Weber <tweber@ill.fr>
 * @date 31-May-2018
 * @license see 'LICENSE' file
 *
 * ----------------------------------------------------------------------------
 * mag-core (part of the Takin software suite)
 * Copyright (C) 2018-2023  Tobias WEBER (Institut Laue-Langevin (ILL),
 *                          Grenoble, France).
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

#ifndef __CLI_WND_H__
#define __CLI_WND_H__


#include <QtCore/QSettings>
#include <QtWidgets/QWidget>
#include <QtWidgets/QDockWidget>
#include <QtWidgets/QComboBox>
#include <QtWidgets/QTextEdit>

#include <memory>

#include "cli/cliparser.h"


/**
 * command line widget
 */
class CommandLineWidget : public QWidget
{
private:
	QSettings *m_pSettings = nullptr;

	QTextEdit *m_pEditHistory = new QTextEdit(this);
	QComboBox *m_pEditCLI = new QComboBox(this);
	QStringList m_completerItems;

	CliParserContext m_parsectx;


protected:
	void CommandEntered();
	void ScrollToEnd();

	void UpdateCompleter();
	void CompleterActivated(const QString &str);

	virtual void resizeEvent(QResizeEvent *evt) override;


public:
	CommandLineWidget(QWidget *pParent = nullptr, QSettings *pSettings = nullptr);
	virtual ~CommandLineWidget();

	CliParserContext& GetParserContext() { return m_parsectx; }

	void PrintOutputString(bool is_err, const QString &str);

	template<typename ...T> void PrintOutput(bool is_err, T&&... msgs)
	{
		std::ostringstream ostr;
		(ostr << ... << std::forward<T>(msgs));
		PrintOutputString(is_err, ostr.str().c_str());
	}

	void SetCompleterItems(const QStringList& lst) { m_completerItems = lst; UpdateCompleter(); }
};



/**
 * the dock which contains the command line widget
 */
class CommandLine : public QDockWidget
{
private:
	std::unique_ptr<CommandLineWidget> m_pCLI;

public:
	CommandLine(QWidget* pParent = nullptr, QSettings *pSettings = nullptr);
	virtual ~CommandLine();

	const CommandLineWidget* GetWidget() const { return m_pCLI.get(); }
	CommandLineWidget* GetWidget() { return m_pCLI.get(); }
};


#endif
