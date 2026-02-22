// This is the SIP interface definition for the QList based mapped types
// specific to the QtGui module.
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


%MappedType QList<QFontDatabase::WritingSystem>
        /TypeHintIn="Sequence[QFontDatabase.WritingSystem]",
        TypeHintOut="List[QFontDatabase.WritingSystem]",
        TypeHintValue="[]"/
{
%TypeHeaderCode
#include <qfontdatabase.h>
%End

%ConvertFromTypeCode
    PyObject *l = PyList_New(sipCpp->size());

    if (!l)
        return 0;

    for (int i = 0; i < sipCpp->size(); ++i)
    {
        PyObject *eobj = sipConvertFromEnum(sipCpp->at(i),
                sipType_QFontDatabase_WritingSystem);

        if (!eobj)
        {
            Py_DECREF(l);

            return 0;
        }

        PyList_SetItem(l, i, eobj);
    }

    return l;
%End

%ConvertToTypeCode
    if (!sipIsErr)
        return (PySequence_Check(sipPy)
#if PY_MAJOR_VERSION < 3
                && !PyString_Check(sipPy)
#endif
                && !PyUnicode_Check(sipPy));

    Py_ssize_t len = PySequence_Size(sipPy);

    if (len < 0)
        return 0;

    QList<QFontDatabase::WritingSystem> *ql = new QList<QFontDatabase::WritingSystem>;

    for (Py_ssize_t i = 0; i < len; ++i)
    {
        PyObject *itm = PySequence_GetItem(sipPy, i);

        if (!itm)
        {
            delete ql;
            *sipIsErr = 1;

            return 0;
        }

        int v = sipConvertToEnum(itm, sipType_QFontDatabase_WritingSystem);

        if (PyErr_Occurred())
        {
            PyErr_Format(PyExc_TypeError,
                    "element %zd has type '%s' but 'QFontDatabase.WritingSystem' is expected",
                    i, sipPyTypeName(Py_TYPE(itm)));

            Py_DECREF(itm);
            delete ql;
            *sipIsErr = 1;

            return 0;
        }

        ql->append(static_cast<QFontDatabase::WritingSystem>(v));

        Py_DECREF(itm);
    }

    *sipCppPtr = ql;

    return sipGetState(sipTransferObj);
%End
};
