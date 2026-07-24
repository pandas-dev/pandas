// This is the SIP interface definition for QAxWidget.
//
// Copyright (c) 2024 Riverbank Computing Limited <info@riverbankcomputing.com>
// 
// This file is part of PyQt5.
// 
// This file may be used under the terms of the GNU General Public License
// version 3.0 as published by the Free Software Foundation and appearing in
// the file LICENSE included in the packaging of this file.  Please review the
// following information to ensure the GNU General Public License version 3.0
// requirements will be met: http://www.gnu.org/copyleft/gpl.html.
// 
// If you do not wish to use this file under the terms of the GPL version 3.0
// then you may purchase a commercial license.  For more information contact
// info@riverbankcomputing.com.
// 
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.


class QAxWidget : QWidget, QAxBase /PyQtNoQMetaObject/
{
%TypeHeaderCode
#include <qaxwidget.h>
%End

public:
    QAxWidget(QWidget *parent /TransferThis/ = 0, Qt::WindowFlags flags = 0);
    QAxWidget(const QString &, QWidget *parent /TransferThis/ = 0,
            Qt::WindowFlags flags = 0);
    //QAxWidget(IUnknown *, QWidget *parent /TransferThis/ = 0, Qt::WindowFlags flags = 0);
    ~QAxWidget();

    void clear();
    bool doVerb(const QString &);

    QSize sizeHint() const;
    QSize minimumSizeHint() const;

    //virtual QaxAggregated *createAggregate();

protected:
    //bool initialize(IUnknown **);
    virtual bool createHostWindow(bool);

    void changeEvent(QEvent *);
    void resizeEvent(QResizeEvent *);

    virtual bool translateKeyEvent(int,int) const;

    void connectNotify(const QMetaMethod &);
};
