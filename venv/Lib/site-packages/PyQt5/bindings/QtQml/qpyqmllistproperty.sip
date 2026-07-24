// This is the SIP specification of the QQmlListProperty mapped type.
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


%MappedType QQmlListProperty<QObject> /TypeHint="QQmlListProperty"/
{
%TypeHeaderCode
#include "qpyqmllistpropertywrapper.h"
%End

%ConvertFromTypeCode
    return qpyqml_QQmlListPropertyWrapper_New(sipCpp, 0);
%End

%ConvertToTypeCode
    if (sipIsErr == NULL)
        return PyObject_IsInstance(sipPy, (PyObject *)qpyqml_QQmlListPropertyWrapper_TypeObject);

    *sipCppPtr = ((qpyqml_QQmlListPropertyWrapper *)sipPy)->qml_list_property;

    // It isn't a temporary copy.
    return 0;
%End
};
