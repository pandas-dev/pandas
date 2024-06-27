// This is the SIP specification of the QPyDBusPendingReply class.
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


class QPyDBusPendingReply : QDBusPendingCall /PyName=QDBusPendingReply/
{
%TypeHeaderCode
#include <qpydbuspendingreply.h>
%End

public:
    QPyDBusPendingReply();
    QPyDBusPendingReply(const QPyDBusPendingReply &other);
    QPyDBusPendingReply(const QDBusPendingCall &call);
    QPyDBusPendingReply(const QDBusMessage &reply);

    // The /ReleaseGIL/ annotation is needed because QDBusPendingCall has an
    // internal mutex.
    QVariant argumentAt(int index) const /ReleaseGIL/;
    QDBusError error() const /ReleaseGIL/;
    bool isError() const /ReleaseGIL/;
    bool isFinished() const /ReleaseGIL/;
    bool isValid() const /ReleaseGIL/;
    QDBusMessage reply() const /ReleaseGIL/;
    void waitForFinished() /ReleaseGIL/;

    SIP_PYOBJECT value(SIP_PYOBJECT type /TypeHintValue="None"/ = 0) const /HoldGIL/;
};
