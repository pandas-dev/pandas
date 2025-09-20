// This is the SIP specification of the qmlAttachedPropertiesObject() function.
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
#include <qqml.h>
%End


QObject *qmlAttachedPropertiesObject(SIP_PYTYPE, QObject *object,
        bool create = true);
%MethodCode
    QObject *proxy = qpyqml_find_proxy_for(a1);

    if (!proxy)
    {
        sipError = sipErrorFail;
    }
    else
    {
        static QHash<PyTypeObject *, int> cache;

        int idx = cache.value((PyTypeObject *)a0, -1);
        const QMetaObject *mo = pyqt5_qtqml_get_qmetaobject((PyTypeObject *)a0);

        sipRes = qmlAttachedPropertiesObject(&idx, proxy, mo, a2);

        cache.insert((PyTypeObject *)a0, idx);
    }
%End
