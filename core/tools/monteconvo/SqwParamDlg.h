/**
 * S(q,w) parameters dialog
 * @author Tobias Weber <tobias.weber@tum.de>
 * @date aug-2015
 * @license GPLv2
 */

#ifndef __SQW_DLG_H__
#define __SQW_DLG_H__

#include <QDialog>
#include <QSettings>

#include "ui/ui_sqwparams.h"
#include "sqw.h"


enum
{
	SQW_NAME = 0,
	SQW_TYPE = 1,
	SQW_VAL = 2,

	SQW_ERR = 3,
	SQW_RANGE = 4,
	SQW_FIT = 5,
};


class SqwParamDlg : public QDialog, Ui::SqwParamDlg
{ Q_OBJECT
protected:
	QSettings *m_pSett = nullptr;

protected:
	void SaveSqwParams();
	virtual void showEvent(QShowEvent *pEvt) override;

public:
	SqwParamDlg(QWidget* pParent=nullptr, QSettings* pSett=nullptr);
	virtual ~SqwParamDlg();

public slots:
	void SqwLoaded(const std::vector<SqwBase::t_var>&, const std::vector<SqwBase::t_var_fit>*);

protected slots:
	void ButtonBoxClicked(QAbstractButton *pBtn);

signals:
	void SqwParamsChanged(const std::vector<SqwBase::t_var>&, const std::vector<SqwBase::t_var_fit>*);
};

#endif
