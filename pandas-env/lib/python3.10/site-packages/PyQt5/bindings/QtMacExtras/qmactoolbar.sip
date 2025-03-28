// This is the SIP interface definition for QMacToolBar.
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


%If (Qt_5_3_0 -)

%If (PyQt_MacOSXOnly)

class NSToolbar;

class QMacToolBar : QObject
{
%TypeHeaderCode
#include <qmactoolbar.h>
%End

public:
    explicit QMacToolBar(QObject *parent /TransferThis/ = 0);
    QMacToolBar(const QString &identifier, QObject *parent /TransferThis/ = 0);
    ~QMacToolBar();

    QMacToolBarItem *addItem(const QIcon &icon, const QString &text);
    QMacToolBarItem *addAllowedItem(const QIcon &icon, const QString &text);
    void addSeparator();

    void setItems(QList<QMacToolBarItem *> &items);
    QList<QMacToolBarItem *> items();
    void setAllowedItems(QList<QMacToolBarItem *> &allowedItems);
    QList<QMacToolBarItem *> allowedItems();

    void attachToWindow(QWindow *window);
    void detachFromWindow();

    NSToolbar *nativeToolbar() const;

private:
    QMacToolBar(const QMacToolBar &);

%ConvertToSubClassCode
    if (sipCpp->inherits(sipName_QMacToolBar))
        sipType = sipType_QMacToolBar;
    else if (sipCpp->inherits(sipName_QMacToolBarItem))
        sipType = sipType_QMacToolBarItem;
    else
        sipType = 0;
%End
};

%End

%End
