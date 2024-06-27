// This is the SIP interface definition for the QVector based mapped types
// specific to the QtGui module.
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


%MappedType QVector<unsigned>
        /TypeHintIn="Iterable[int]", TypeHintOut="List[int]",
        TypeHintValue="[]"/
{
%TypeHeaderCode
#include <qvector.h>
%End

%ConvertFromTypeCode
    PyObject *l = PyList_New(sipCpp->size());

    if (!l)
        return 0;

    for (int i = 0; i < sipCpp->size(); ++i)
    {
        // Convert to a Python long to make sure it doesn't get interpreted as
        // a signed value.
        PyObject *pobj = PyLong_FromUnsignedLong(sipCpp->value(i));

        if (!pobj)
        {
            Py_DECREF(l);

            return 0;
        }

        PyList_SetItem(l, i, pobj);
    }

    return l;
%End

%ConvertToTypeCode
    PyObject *iter = PyObject_GetIter(sipPy);

    if (!sipIsErr)
    {
        PyErr_Clear();
        Py_XDECREF(iter);

        return (iter
#if PY_MAJOR_VERSION < 3
                && !PyString_Check(sipPy)
#endif
                && !PyUnicode_Check(sipPy));
    }

    if (!iter)
    {
        *sipIsErr = 1;

        return 0;
    }

    QVector<unsigned> *qv = new QVector<unsigned>;

    for (Py_ssize_t i = 0; ; ++i)
    {
        PyErr_Clear();
        PyObject *itm = PyIter_Next(iter);

        if (!itm)
        {
            if (PyErr_Occurred())
            {
                delete qv;
                Py_DECREF(iter);
                *sipIsErr = 1;

                return 0;
            }

            break;
        }

        PyErr_Clear();
        unsigned long val = PyLong_AsUnsignedLongMask(itm);

        if (PyErr_Occurred())
        {
            PyErr_Format(PyExc_TypeError,
                    "index %zd has type '%s' but 'int' is expected", i,
                    sipPyTypeName(Py_TYPE(itm)));

            Py_DECREF(itm);
            delete qv;
            Py_DECREF(iter);
            *sipIsErr = 1;

            return 0;
        }

        qv->append(val);

        Py_DECREF(itm);
    }

    Py_DECREF(iter);

    *sipCppPtr = qv;

    return sipGetState(sipTransferObj);
%End
};


%MappedType QVector<float>
        /TypeHintIn="Iterable[float]", TypeHintOut="List[float]",
        TypeHintValue="[]"/
{
%TypeHeaderCode
#include <qvector.h>
%End

%ConvertFromTypeCode
    PyObject *l = PyList_New(sipCpp->size());

    if (!l)
        return 0;

    for (int i = 0; i < sipCpp->size(); ++i)
    {
        PyObject *pobj = PyFloat_FromDouble(sipCpp->value(i));

        if (!pobj)
        {
            Py_DECREF(l);

            return 0;
        }

        PyList_SetItem(l, i, pobj);
    }

    return l;
%End

%ConvertToTypeCode
    PyObject *iter = PyObject_GetIter(sipPy);

    if (!sipIsErr)
    {
        PyErr_Clear();
        Py_XDECREF(iter);

        return (iter
#if PY_MAJOR_VERSION < 3
                && !PyString_Check(sipPy)
#endif
                && !PyUnicode_Check(sipPy));
    }

    if (!iter)
    {
        *sipIsErr = 1;

        return 0;
    }

    QVector<float> *qv = new QVector<float>;

    for (Py_ssize_t i = 0; ; ++i)
    {
        PyErr_Clear();
        PyObject *itm = PyIter_Next(iter);

        if (!itm)
        {
            if (PyErr_Occurred())
            {
                delete qv;
                Py_DECREF(iter);
                *sipIsErr = 1;

                return 0;
            }

            break;
        }

        PyErr_Clear();
        double val = PyFloat_AsDouble(itm);

        if (PyErr_Occurred())
        {
            PyErr_Format(PyExc_TypeError,
                    "index %zd has type '%s' but 'int' is expected", i,
                    sipPyTypeName(Py_TYPE(itm)));

            Py_DECREF(itm);
            delete qv;
            Py_DECREF(iter);
            *sipIsErr = 1;

            return 0;
        }

        qv->append(val);

        Py_DECREF(itm);
    }

    Py_DECREF(iter);

    *sipCppPtr = qv;

    return sipGetState(sipTransferObj);
%End
};


%If (PyQt_qreal_double)

%MappedType QVector<qreal>
        /TypeHintIn="Iterable[float]", TypeHintOut="List[float]",
        TypeHintValue="[]"/
{
%TypeHeaderCode
#include <qvector.h>
%End

%ConvertFromTypeCode
    PyObject *l = PyList_New(sipCpp->size());

    if (!l)
        return 0;

    for (int i = 0; i < sipCpp->size(); ++i)
    {
        PyObject *pobj = PyFloat_FromDouble(sipCpp->value(i));

        if (!pobj)
        {
            Py_DECREF(l);

            return 0;
        }

        PyList_SetItem(l, i, pobj);
    }

    return l;
%End

%ConvertToTypeCode
    PyObject *iter = PyObject_GetIter(sipPy);

    if (!sipIsErr)
    {
        PyErr_Clear();
        Py_XDECREF(iter);

        return (iter
#if PY_MAJOR_VERSION < 3
                && !PyString_Check(sipPy)
#endif
                && !PyUnicode_Check(sipPy));
    }

    if (!iter)
    {
        *sipIsErr = 1;

        return 0;
    }

    QVector<qreal> *qv = new QVector<qreal>;

    for (Py_ssize_t i = 0; ; ++i)
    {
        PyErr_Clear();
        PyObject *itm = PyIter_Next(iter);

        if (!itm)
        {
            if (PyErr_Occurred())
            {
                delete qv;
                Py_DECREF(iter);
                *sipIsErr = 1;

                return 0;
            }

            break;
        }

        PyErr_Clear();
        double val = PyFloat_AsDouble(itm);

        if (PyErr_Occurred())
        {
            PyErr_Format(PyExc_TypeError,
                    "index %zd has type '%s' but 'int' is expected", i,
                    sipPyTypeName(Py_TYPE(itm)));

            Py_DECREF(itm);
            delete qv;
            Py_DECREF(iter);
            *sipIsErr = 1;

            return 0;
        }

        qv->append(val);

        Py_DECREF(itm);
    }

    Py_DECREF(iter);

    *sipCppPtr = qv;

    return sipGetState(sipTransferObj);
%End
};

%End


%If (PyQt_Desktop_OpenGL)

%MappedType QVector<GLuint64>
        /TypeHintIn="Iterable[int]", TypeHintOut="List[int]",
        TypeHintValue="[]"/
{
%TypeHeaderCode
#include <qvector.h>
#include <qopengl.h>
%End

%ConvertFromTypeCode
    PyObject *l = PyList_New(sipCpp->size());

    if (!l)
        return 0;

    for (int i = 0; i < sipCpp->size(); ++i)
    {
        // Convert to a Python long to make sure it doesn't get interpreted as
        // a signed value.
        PyObject *pobj = PyLong_FromUnsignedLongLong(sipCpp->value(i));

        if (!pobj)
        {
            Py_DECREF(l);

            return 0;
        }

        PyList_SetItem(l, i, pobj);
    }

    return l;
%End

%ConvertToTypeCode
    PyObject *iter = PyObject_GetIter(sipPy);

    if (!sipIsErr)
    {
        PyErr_Clear();
        Py_XDECREF(iter);

        return (iter
#if PY_MAJOR_VERSION < 3
                && !PyString_Check(sipPy)
#endif
                && !PyUnicode_Check(sipPy));
    }

    if (!iter)
    {
        *sipIsErr = 1;

        return 0;
    }

    QVector<GLuint64> *qv = new QVector<GLuint64>;

    for (Py_ssize_t i = 0; ; ++i)
    {
        PyErr_Clear();
        PyObject *itm = PyIter_Next(iter);

        if (!itm)
        {
            if (PyErr_Occurred())
            {
                delete qv;
                Py_DECREF(iter);
                *sipIsErr = 1;

                return 0;
            }

            break;
        }

        PyErr_Clear();
        unsigned PY_LONG_LONG val = PyLong_AsUnsignedLongLongMask(itm);

        if (PyErr_Occurred())
        {
            PyErr_Format(PyExc_TypeError,
                    "index %zd has type '%s' but 'int' is expected", i,
                    sipPyTypeName(Py_TYPE(itm)));

            Py_DECREF(itm);
            delete qv;
            Py_DECREF(iter);
            *sipIsErr = 1;

            return 0;
        }

        qv->append(val);

        Py_DECREF(itm);
    }

    Py_DECREF(iter);

    *sipCppPtr = qv;

    return sipGetState(sipTransferObj);
%End
};

%End
