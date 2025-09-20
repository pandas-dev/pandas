// This is the SIP interface definition for the quint128 mapped type.
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


%MappedType quint128 /TypeHint="Tuple[int, int, int, int, int, int, int, int, int, int, int, int, int, int, int, int]"/
{
%TypeHeaderCode
#include <qbluetoothuuid.h>
%End

%ConvertFromTypeCode
    PyObject *t = PyTuple_New(16);

    if (!t)
        return 0;

    for (Py_ssize_t i = 0; i < 16; ++i)
    {
        // Convert to a Python long to make sure it doesn't get interpreted as
        // a signed value.
        PyObject *pobj = PyLong_FromUnsignedLong(sipCpp->data[i]);

        if (!pobj)
        {
            Py_DECREF(t);

            return 0;
        }

        PyTuple_SetItem(t, i, pobj);
    }

    return t;
%End

%ConvertToTypeCode
    if (!sipIsErr)
        return (PySequence_Check(sipPy)
#if PY_MAJOR_VERSION < 3
                && !PyString_Check(sipPy)
#endif
                && !PyUnicode_Check(sipPy));

    Py_ssize_t len = PySequence_Size(sipPy);

    if (len != 16)
    {
        // A negative length should only be an internal error so let the
        // original exception stand.
        if (len >= 0)
            PyErr_Format(PyExc_TypeError,
                    "sequence has %zd elements but 16 elements are expected",
                    len);

        *sipIsErr = 1;

        return 0;
    }

    quint128 *qv = new quint128;

    for (Py_ssize_t i = 0; i < 16; ++i)
    {
        PyObject *itm = PySequence_GetItem(sipPy, i);

        if (!itm)
        {
            delete qv;
            *sipIsErr = 1;

            return 0;
        }

        PyErr_Clear();
        unsigned long val = PyLong_AsUnsignedLongMask(itm);

        if (PyErr_Occurred())
        {
            PyErr_Format(PyExc_TypeError,
                    "element %zd has type '%s' but 'int' is expected", i,
                    sipPyTypeName(Py_TYPE(itm)));

            Py_DECREF(itm);
            delete qv;
            *sipIsErr = 1;

            return 0;
        }

        qv->data[i] = val;

        Py_DECREF(itm);
    }

    *sipCppPtr = qv;

    return sipGetState(sipTransferObj);
%End
};
