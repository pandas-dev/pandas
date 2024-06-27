// This is the SIP interface definition for the majority of the QSet based
// mapped types.
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


template<_TYPE_>
%MappedType QSet<_TYPE_>
        /TypeHintIn="Iterable[_TYPE_]", TypeHintOut="Set[_TYPE_]",
        TypeHintValue="set()"/
{
%TypeHeaderCode
#include <qset.h>
%End

%ConvertFromTypeCode
    PyObject *s = PySet_New(0);

    if (!s)
        return 0;

    QSet<_TYPE_>::const_iterator it = sipCpp->constBegin();
    QSet<_TYPE_>::const_iterator end = sipCpp->constEnd();

    while (it != end)
    {
        _TYPE_ *t = new _TYPE_(*it);
        PyObject *tobj = sipConvertFromNewType(t, sipType__TYPE_,
                sipTransferObj);

        if (!tobj)
        {
            delete t;
            Py_DECREF(s);

            return 0;
        }

        PySet_Add(s, tobj);

        ++it;
    }

    return s;
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

    QSet<_TYPE_> *qs = new QSet<_TYPE_>;
 
    for (Py_ssize_t i = 0; ; ++i)
    {
        PyErr_Clear();
        PyObject *itm = PyIter_Next(iter);

        if (!itm)
        {
            if (PyErr_Occurred())
            {
                delete qs;
                Py_DECREF(iter);
                *sipIsErr = 1;

                return 0;
            }

            break;
        }

        int state;
        _TYPE_ *t = reinterpret_cast<_TYPE_ *>(
                sipForceConvertToType(itm, sipType__TYPE_, sipTransferObj,
                        SIP_NOT_NONE, &state, sipIsErr));

        if (*sipIsErr)
        {
            PyErr_Format(PyExc_TypeError,
                    "index %zd has type '%s' but '_TYPE_' is expected", i,
                    sipPyTypeName(Py_TYPE(itm)));

            Py_DECREF(itm);
            delete qs;
            Py_DECREF(iter);

            return 0;
        }

        qs->insert(*t);

        sipReleaseType(t, sipType__TYPE_, state);
        Py_DECREF(itm);
    }
 
    Py_DECREF(iter);

    *sipCppPtr = qs;
 
    return sipGetState(sipTransferObj);
%End
};


template<_TYPE_>
%MappedType QSet<_TYPE_ *>
        /TypeHintIn="Iterable[_TYPE_]", TypeHintOut="Set[_TYPE_]",
        TypeHintValue="set()"/
{
%TypeHeaderCode
#include <qset.h>
%End

%ConvertFromTypeCode
    int gc_enabled = sipEnableGC(0);
    PyObject *s = PySet_New(0);

    if (s)
    {
        QSet<_TYPE_ *>::const_iterator it = sipCpp->constBegin();
        QSet<_TYPE_ *>::const_iterator end = sipCpp->constEnd();

        while (it != end)
        {
            // The explicit (void *) cast allows _TYPE_ to be const.
            PyObject *tobj = sipConvertFromType((void *)*it, sipType__TYPE_,
                    sipTransferObj);

            if (!tobj)
            {
                Py_DECREF(s);
                s = 0;

                break;
            }

            PySet_Add(s, tobj);

            ++it;
        }
    }

    sipEnableGC(gc_enabled);

    return s;
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

    QSet<_TYPE_ *> *qs = new QSet<_TYPE_ *>;
 
    for (Py_ssize_t i = 0; ; ++i)
    {
        PyErr_Clear();
        PyObject *itm = PyIter_Next(iter);

        if (!itm)
        {
            if (PyErr_Occurred())
            {
                delete qs;
                Py_DECREF(iter);
                *sipIsErr = 1;

                return 0;
            }

            break;
        }

        _TYPE_ *t = reinterpret_cast<_TYPE_ *>(
                sipForceConvertToType(itm, sipType__TYPE_, sipTransferObj, 0,
                        0, sipIsErr));

        if (*sipIsErr)
        {
            PyErr_Format(PyExc_TypeError,
                    "index %zd has type '%s' but '_TYPE_' is expected", i,
                    sipPyTypeName(Py_TYPE(itm)));

            Py_DECREF(itm);
            delete qs;
            Py_DECREF(iter);

            return 0;
        }

        qs->insert(t);

        Py_DECREF(itm);
    }
 
    Py_DECREF(iter);

    *sipCppPtr = qs;
 
    return sipGetState(sipTransferObj);
%End
};
