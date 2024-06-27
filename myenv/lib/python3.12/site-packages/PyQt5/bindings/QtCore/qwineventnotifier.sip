// This is the SIP specification of the QWinEventNotifier class.
//
// Copyright (c) 2023 Riverbank Computing Limited <info@riverbankcomputing.com>
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


%If (WS_WIN)

// This hack is for the activated() signal.
typedef Qt::HANDLE HANDLE;

class QWinEventNotifier: QObject
{
%TypeHeaderCode
#include <qwineventnotifier.h>
%End

public:
    explicit QWinEventNotifier(QObject *parent /TransferThis/ = 0);
    explicit QWinEventNotifier(Qt::HANDLE hEvent, QObject *parent /TransferThis/ = 0);
    ~QWinEventNotifier();

    Qt::HANDLE handle() const;
    bool isEnabled() const;
    void setHandle(Qt::HANDLE hEvent);

public slots:
    void setEnabled(bool enable);

signals:
    void activated(HANDLE hEvent);

protected:
    bool event(QEvent *e);

private:
    QWinEventNotifier(const QWinEventNotifier &);
};

%End
