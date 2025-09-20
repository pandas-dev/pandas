// This is the SIP specification of the qmlRegisterType() function.
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


%ModuleHeaderCode
#include <qpyqml_api.h>
%End


%ModuleCode
// Imports from QtCore.
pyqt5_qtqml_get_qmetaobject_t pyqt5_qtqml_get_qmetaobject;
%End


%PostInitialisationCode
// Imports from QtCore.
pyqt5_qtqml_get_qmetaobject = (pyqt5_qtqml_get_qmetaobject_t)sipImportSymbol(
        "pyqt5_get_qmetaobject");
Q_ASSERT(pyqt5_qtqml_get_qmetaobject);
%End


int qmlRegisterRevision(SIP_PYTYPE, int revision, const char *uri, int major,
        int minor, SIP_PYTYPE attachedProperties = 0);
%MethodCode
    if ((sipRes = qpyqml_register_library_type((PyTypeObject *)a0, a2, a3, a4, 0, a1, (PyTypeObject *)a5)) < 0)
        sipError = sipErrorFail;
%End


%If (Qt_5_2_0 -)
int qmlRegisterSingletonType(const QUrl &url, const char *uri, int major,
        int minor, const char *qmlName);
%End


int qmlRegisterSingletonType(SIP_PYTYPE, const char *uri, int major, int minor,
        const char *typeName, SIP_PYCALLABLE factory /TypeHint="Callable[[QQmlEngine, QJSEngine], Any]"/);
%MethodCode
    if ((sipRes = qpyqml_register_singleton_type((PyTypeObject *)a0, a1, a2, a3, a4, a5)) < 0)
        sipError = sipErrorFail;
%End


int qmlRegisterType(const QUrl &url, const char *uri, int major, int minor,
        const char *qmlName);


int qmlRegisterType(SIP_PYTYPE, SIP_PYTYPE attachedProperties = 0);
%MethodCode
    if ((sipRes = qpyqml_register_type((PyTypeObject *)a0, (PyTypeObject *)a1)) < 0)
        sipError = sipErrorFail;
%End


int qmlRegisterType(SIP_PYTYPE, const char *uri, int major, int minor,
        const char *qmlName, SIP_PYTYPE attachedProperties = 0);
%MethodCode
    if ((sipRes = qpyqml_register_library_type((PyTypeObject *)a0, a1, a2, a3, a4, -1, (PyTypeObject *)a5)) < 0)
        sipError = sipErrorFail;
%End


int qmlRegisterType(SIP_PYTYPE, int revision, const char *uri, int major,
        int minor, const char *qmlName, SIP_PYTYPE attachedProperties = 0);
%MethodCode
    if ((sipRes = qpyqml_register_library_type((PyTypeObject *)a0, a2, a3, a4, a5, a1, (PyTypeObject *)a6)) < 0)
        sipError = sipErrorFail;
%End


int qmlRegisterUncreatableType(SIP_PYTYPE, const char *uri, int major,
        int minor, const char *qmlName, const QString &reason);
%MethodCode
    if ((sipRes = qpyqml_register_uncreatable_type((PyTypeObject *)a0, a1, a2, a3, a4, *a5, -1)) < 0)
        sipError = sipErrorFail;
%End


int qmlRegisterUncreatableType(SIP_PYTYPE, int revision, const char *uri,
        int major, int minor, const char *qmlName, const QString &reason);
%MethodCode
    if ((sipRes = qpyqml_register_uncreatable_type((PyTypeObject *)a0, a2, a3, a4, a5, *a6, a1)) < 0)
        sipError = sipErrorFail;
%End
