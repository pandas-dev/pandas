// This is the SIP interface definition for the majority of the QVector based
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
%MappedType QVector<_TYPE_>
        /TypeHintIn="Iterable[_TYPE_]", TypeHintOut="List[_TYPE_]",
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
        _TYPE_ *t = new _TYPE_(sipCpp->at(i));
        PyObject *tobj = sipConvertFromNewType(t, sipType__TYPE_,
                sipTransferObj);

        if (!tobj)
        {
            delete t;
            Py_DECREF(l);

            return 0;
        }

        PyList_SetItem(l, i, tobj);
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

    QVector<_TYPE_> *qv = new QVector<_TYPE_>;
 
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
            delete qv;
            Py_DECREF(iter);

            return 0;
        }

        qv->append(*t);

        sipReleaseType(t, sipType__TYPE_, state);
        Py_DECREF(itm);
    }

    Py_DECREF(iter);
 
    *sipCppPtr = qv;
 
    return sipGetState(sipTransferObj);
%End
};


template<_TYPE_>
%MappedType QVector<_TYPE_ *>
        /TypeHintIn="Iterable[_TYPE_]", TypeHintOut="List[_TYPE_]",
        TypeHintValue="[]"/
{
%TypeHeaderCode
#include <qvector.h>
%End

%ConvertFromTypeCode
    int gc_enabled = sipEnableGC(0);
    PyObject *l = PyList_New(sipCpp->size());

    if (l)
    {
        for (int i = 0; i < sipCpp->size(); ++i)
        {
            _TYPE_ *t = sipCpp->at(i);

            // The explicit (void *) cast allows _TYPE_ to be const.
            PyObject *tobj = sipConvertFromNewType((void *)t, sipType__TYPE_,
                    sipTransferObj);

            if (!tobj)
            {
                Py_DECREF(l);
                l = 0;

                break;
            }

            PyList_SetItem(l, i, tobj);
        }
    }

    sipEnableGC(gc_enabled);

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

    QVector<_TYPE_ *> *qv = new QVector<_TYPE_ *>;
 
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

        _TYPE_ *t = reinterpret_cast<_TYPE_ *>(
                sipForceConvertToType(itm, sipType__TYPE_, sipTransferObj, 0,
                        0, sipIsErr));

        if (*sipIsErr)
        {
            PyErr_Format(PyExc_TypeError,
                    "index %zd has type '%s' but '_TYPE_' is expected", i,
                    sipPyTypeName(Py_TYPE(itm)));

            Py_DECREF(itm);
            delete qv;
            Py_DECREF(iter);

            return 0;
        }

        qv->append(t);

        Py_DECREF(itm);
    }

    Py_DECREF(iter);
 
    *sipCppPtr = qv;
 
    return sipGetState(sipTransferObj);
%End
};


template<qreal, _TYPE_>
%MappedType QVector<QPair<qreal, _TYPE_> >
        /TypeHintIn="Iterable[Tuple[float, _TYPE_]]",
        TypeHintOut="List[Tuple[float, _TYPE_]]", TypeHintValue="[]"/
{
%TypeHeaderCode
#include <qvector.h>
#include <qpair.h>
%End

%ConvertFromTypeCode
    PyObject *l = PyList_New(sipCpp->size());

    if (!l)
        return 0;

    for (int i = 0; i < sipCpp->size(); ++i)
    {
        const QPair<qreal, _TYPE_> &p = sipCpp->at(i);
        _TYPE_ *s2 = new _TYPE_(p.second);
        PyObject *pobj = sipBuildResult(NULL, "(dN)", (double)p.first, s2,
                sipType__TYPE_, sipTransferObj);

        if (!pobj)
        {
            delete s2;
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

    QVector<QPair<qreal, _TYPE_> > *qv = new QVector<QPair<qreal, _TYPE_> >;
 
    for (Py_ssize_t i = 0; ; ++i)
    {
        PyErr_Clear();
        PyObject *seq = PyIter_Next(iter);

        if (!seq)
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

        Py_ssize_t sub_len;

        if (PySequence_Check(seq)
#if PY_MAJOR_VERSION < 3
                && !PyString_Check(seq)
#endif
                && !PyUnicode_Check(seq))
            sub_len = PySequence_Size(seq);
        else
            sub_len = -1;

        if (sub_len != 2)
        {
            if (sub_len < 0)
                PyErr_Format(PyExc_TypeError,
                        "index %zd has type '%s' but a 2 element non-string sequence is expected",
                        i, sipPyTypeName(Py_TYPE(seq)));
            else
                PyErr_Format(PyExc_TypeError,
                        "index %zd is a sequence of %zd sub-elements but 2 sub-elements are expected",
                        i, sub_len);

            Py_DECREF(seq);
            delete qv;
            Py_DECREF(iter);
            *sipIsErr = 1;

            return 0;
        }

        PyObject *itm1 = PySequence_GetItem(seq, 0);

        if (!itm1)
        {
            Py_DECREF(seq);
            delete qv;
            Py_DECREF(iter);
            *sipIsErr = 1;

            return 0;
        }

        PyErr_Clear();
        qreal s1 = PyFloat_AsDouble(itm1);

        if (PyErr_Occurred())
        {
            PyErr_Format(PyExc_TypeError,
                    "the first sub-element of index %zd has type '%s' but 'float' is expected",
                    i, sipPyTypeName(Py_TYPE(itm1)));

            Py_DECREF(itm1);
            Py_DECREF(seq);
            delete qv;
            Py_DECREF(iter);
            *sipIsErr = 1;

            return 0;
        }

        PyObject *itm2 = PySequence_GetItem(seq, 1);

        if (!itm2)
        {
            Py_DECREF(itm1);
            Py_DECREF(seq);
            delete qv;
            Py_DECREF(iter);
            *sipIsErr = 1;

            return 0;
        }

        int state2;
        _TYPE_ *s2 = reinterpret_cast<_TYPE_ *>(
                sipForceConvertToType(itm2, sipType__TYPE_, sipTransferObj,
                        SIP_NOT_NONE, &state2, sipIsErr));
 
        if (*sipIsErr)
        {
            PyErr_Format(PyExc_TypeError,
                    "the second sub-element of index %zd has type '%s' but '_TYPE_' is expected",
                    i, sipPyTypeName(Py_TYPE(itm2)));

            Py_DECREF(itm2);
            Py_DECREF(itm1);
            Py_DECREF(seq);
            delete qv;
            Py_DECREF(iter);

            return 0;
        }

        qv->append(QPair<qreal, _TYPE_>(s1, *s2));

        sipReleaseType(s2, sipType__TYPE_, state2);
        Py_DECREF(itm2);
        Py_DECREF(itm1);
        Py_DECREF(seq);
    }

    Py_DECREF(iter);
 
    *sipCppPtr = qv;
 
    return sipGetState(sipTransferObj);
%End
};


%MappedType QVector<int>
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
        PyObject *pobj = SIPLong_FromLong(sipCpp->value(i));

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

    QVector<int> *qv = new QVector<int>;
 
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

        int val = sipLong_AsInt(itm);

        if (PyErr_Occurred())
        {
            if (PyErr_ExceptionMatches(PyExc_TypeError))
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


%If (Qt_5_12_0 -)

%MappedType QVector<quint16>
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
        PyObject *pobj = SIPLong_FromLong(sipCpp->value(i));

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

    QVector<quint16> *qv = new QVector<quint16>;
 
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

        quint16 val = sipLong_AsUnsignedShort(itm);

        if (PyErr_Occurred())
        {
            if (PyErr_ExceptionMatches(PyExc_TypeError))
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
