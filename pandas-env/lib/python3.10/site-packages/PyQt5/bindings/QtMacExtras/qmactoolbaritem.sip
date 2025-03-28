// This is the SIP interface definition for QMacToolBarItem.
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

class NSToolbarItem;

class QMacToolBarItem : QObject
{
%TypeHeaderCode
#include <qmactoolbaritem.h>
%End

public:
    enum StandardItem
    {
        NoStandardItem,
        Space,
        FlexibleSpace
    };

    QMacToolBarItem(QObject *parent /TransferThis/ = 0);
    virtual ~QMacToolBarItem();

    bool selectable() const;
    void setSelectable(bool selectable);

    StandardItem standardItem() const;
    void setStandardItem(StandardItem standardItem);

    QString text() const;
    void setText(const QString &text);

    QIcon icon() const;
    void setIcon(const QIcon &icon);

    NSToolbarItem *nativeToolBarItem() const;

signals:
    void activated();

private:
    QMacToolBarItem(const QMacToolBarItem &);
};

%End

%End
