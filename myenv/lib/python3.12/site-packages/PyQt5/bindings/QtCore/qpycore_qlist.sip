// This is the SIP interface definition for the majority of the QList based
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


// Note that when we test the type of an object to see if it can be converted
// to a collection we only check if it is iterable.  We do not check the
// types of the contents - we assume we will be able to convert them when
// requested.  This allows us to raise exceptions specific to an individual
// object.  This approach doesn't work if there are overloads that can only be
// distinguished by the types of the template arguments.  Currently there are
// no such cases in PyQt5.  Note also that we don't consider strings to be
// iterables.


template<_TYPE_>
%MappedType QList<_TYPE_>
        /TypeHintIn="Iterable[_TYPE_]", TypeHintOut="List[_TYPE_]",
        TypeHintValue="[]"/
{
%TypeHeaderCode
#include <qlist.h>
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

    QList<_TYPE_> *ql = new QList<_TYPE_>;
 
    for (Py_ssize_t i = 0; ; ++i)
    {
        PyErr_Clear();
        PyObject *itm = PyIter_Next(iter);

        if (!itm)
        {
            if (PyErr_Occurred())
            {
                delete ql;
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
            delete ql;
            Py_DECREF(iter);

            return 0;
        }

        ql->append(*t);

        sipReleaseType(t, sipType__TYPE_, state);
        Py_DECREF(itm);
    }
 
    Py_DECREF(iter);

    *sipCppPtr = ql;
 
    return sipGetState(sipTransferObj);
%End
};


template<_TYPE_>
%MappedType QList<_TYPE_ *>
        /TypeHintIn="Iterable[_TYPE_]", TypeHintOut="List[_TYPE_]",
        TypeHintValue="[]"/
{
%TypeHeaderCode
#include <qlist.h>
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
            PyObject *tobj = sipConvertFromType((void *)t, sipType__TYPE_,
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

    QList<_TYPE_ *> *ql = new QList<_TYPE_ *>;
 
    for (Py_ssize_t i = 0; ; ++i)
    {
        PyErr_Clear();
        PyObject *itm = PyIter_Next(iter);

        if (!itm)
        {
            if (PyErr_Occurred())
            {
                delete ql;
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
            delete ql;
            Py_DECREF(iter);

            return 0;
        }

        ql->append(t);

        Py_DECREF(itm);
    }

    Py_DECREF(iter);
 
    *sipCppPtr = ql;
 
    return sipGetState(sipTransferObj);
%End
};


template<_TYPE1_, _TYPE2_>
%MappedType QList<QPair<_TYPE1_, _TYPE2_> >
        /TypeHintIn="Iterable[Tuple[_TYPE1_, _TYPE2_]]",
        TypeHintOut="List[Tuple[_TYPE1_, _TYPE2_]]", TypeHintValue="[]"/
{
%TypeHeaderCode
#include <qlist.h>
#include <qpair.h>
%End

%ConvertFromTypeCode
    PyObject *l = PyList_New(sipCpp->size());

    if (!l)
        return 0;

    for (int i = 0; i < sipCpp->size(); ++i)
    {
        const QPair<_TYPE1_, _TYPE2_> &p = sipCpp->at(i);
        _TYPE1_ *s1 = new _TYPE1_(p.first);
        _TYPE2_ *s2 = new _TYPE2_(p.second);
        PyObject *pobj = sipBuildResult(NULL, "(NN)", s1, sipType__TYPE1_,
                sipTransferObj, s2, sipType__TYPE2_, sipTransferObj);

        if (!pobj)
        {
            delete s1;
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

    QList<QPair<_TYPE1_, _TYPE2_> > *ql = new QList<QPair<_TYPE1_, _TYPE2_> >;
 
    for (Py_ssize_t i = 0; ; ++i)
    {
        PyErr_Clear();
        PyObject *seq = PyIter_Next(iter);

        if (!seq)
        {
            if (PyErr_Occurred())
            {
                delete ql;
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
            delete ql;
            Py_DECREF(iter);
            *sipIsErr = 1;

            return 0;
        }

        PyObject *itm1 = PySequence_GetItem(seq, 0);

        if (!itm1)
        {
            Py_DECREF(seq);
            delete ql;
            Py_DECREF(iter);
            *sipIsErr = 1;

            return 0;
        }

        int state1;
        _TYPE1_ *s1 = reinterpret_cast<_TYPE1_ *>(
                sipForceConvertToType(itm1, sipType__TYPE1_, sipTransferObj,
                        SIP_NOT_NONE, &state1, sipIsErr));

        if (*sipIsErr)
        {
            PyErr_Format(PyExc_TypeError,
                    "the first sub-element of index %zd has type '%s' but '_TYPE1_' is expected",
                    i, sipPyTypeName(Py_TYPE(itm1)));

            Py_DECREF(itm1);
            Py_DECREF(seq);
            delete ql;
            Py_DECREF(iter);

            return 0;
        }

        PyObject *itm2 = PySequence_GetItem(seq, 1);

        if (!itm2)
        {
            sipReleaseType(s1, sipType__TYPE1_, state1);
            Py_DECREF(itm1);
            Py_DECREF(seq);
            delete ql;
            Py_DECREF(iter);
            *sipIsErr = 1;

            return 0;
        }

        int state2;
        _TYPE2_ *s2 = reinterpret_cast<_TYPE2_ *>(
                sipForceConvertToType(itm2, sipType__TYPE2_, sipTransferObj,
                        SIP_NOT_NONE, &state2, sipIsErr));
 
        if (*sipIsErr)
        {
            PyErr_Format(PyExc_TypeError,
                    "the second sub-element of index %zd has type '%s' but '_TYPE2_' is expected",
                    i, sipPyTypeName(Py_TYPE(itm2)));

            Py_DECREF(itm2);
            sipReleaseType(s1, sipType__TYPE1_, state1);
            Py_DECREF(itm1);
            Py_DECREF(seq);
            delete ql;
            Py_DECREF(iter);

            return 0;
        }

        ql->append(QPair<_TYPE1_, _TYPE2_>(*s1, *s2));

        sipReleaseType(s2, sipType__TYPE2_, state2);
        Py_DECREF(itm2);
        sipReleaseType(s1, sipType__TYPE1_, state1);
        Py_DECREF(itm1);
        Py_DECREF(seq);
    }

    Py_DECREF(iter);
 
    *sipCppPtr = ql;
 
    return sipGetState(sipTransferObj);
%End
};


%If (Qt_5_1_0 -)

%MappedType QList<QPair<int, int> >
        /TypeHintIn="Iterable[Tuple[int, int]]",
        TypeHintOut="List[Tuple[int, int]]", TypeHintValue="[]"/
{
%TypeHeaderCode
#include <qlist.h>
#include <qpair.h>
%End

%ConvertFromTypeCode
    PyObject *l = PyList_New(sipCpp->size());

    if (!l)
        return 0;

    for (int i = 0; i < sipCpp->size(); ++i)
    {
        const QPair<int, int> &p = sipCpp->at(i);
        PyObject *pobj = Py_BuildValue((char *)"ii", p.first, p.second);

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

    QList<QPair<int, int> > *ql = new QList<QPair<int, int> >;
 
    for (Py_ssize_t i = 0; ; ++i)
    {
        PyErr_Clear();
        PyObject *seq = PyIter_Next(iter);

        if (!seq)
        {
            if (PyErr_Occurred())
            {
                delete ql;
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
            delete ql;
            Py_DECREF(iter);
            *sipIsErr = 1;

            return 0;
        }

        PyObject *itm1 = PySequence_GetItem(seq, 0);

        if (!itm1)
        {
            Py_DECREF(seq);
            delete ql;
            Py_DECREF(iter);
            *sipIsErr = 1;

            return 0;
        }

        int first = sipLong_AsInt(itm1);

        if (PyErr_Occurred())
        {
            if (PyErr_ExceptionMatches(PyExc_TypeError))
                PyErr_Format(PyExc_TypeError,
                        "the first sub-element of index %zd has type '%s' but 'int' is expected",
                        i, sipPyTypeName(Py_TYPE(itm1)));

            Py_DECREF(itm1);
            Py_DECREF(seq);
            delete ql;
            Py_DECREF(iter);
            *sipIsErr = 1;

            return 0;
        }

        PyObject *itm2 = PySequence_GetItem(seq, 1);

        if (!itm2)
        {
            Py_DECREF(itm1);
            Py_DECREF(seq);
            delete ql;
            Py_DECREF(iter);
            *sipIsErr = 1;

            return 0;
        }

        int second = sipLong_AsInt(itm2);

        if (PyErr_Occurred())
        {
            if (PyErr_ExceptionMatches(PyExc_TypeError))
                PyErr_Format(PyExc_TypeError,
                        "the second sub-element of index %zd has type '%s' but 'int' is expected",
                        i, sipPyTypeName(Py_TYPE(itm2)));

            Py_DECREF(itm2);
            Py_DECREF(itm1);
            Py_DECREF(seq);
            delete ql;
            Py_DECREF(iter);
            *sipIsErr = 1;

            return 0;
        }

        ql->append(QPair<int, int>(first, second));

        Py_DECREF(itm2);
        Py_DECREF(itm1);
        Py_DECREF(seq);
    }

    Py_DECREF(iter);
 
    *sipCppPtr = ql;
 
    return sipGetState(sipTransferObj);
%End
};

%End


%MappedType QList<int>
        /TypeHintIn="Iterable[int]", TypeHintOut="List[int]",
        TypeHintValue="[]"/
{
%TypeHeaderCode
#include <qlist.h>
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

    QList<int> *ql = new QList<int>;
 
    for (Py_ssize_t i = 0; ; ++i)
    {
        PyErr_Clear();
        PyObject *itm = PyIter_Next(iter);

        if (!itm)
        {
            if (PyErr_Occurred())
            {
                delete ql;
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
            delete ql;
            Py_DECREF(iter);
            *sipIsErr = 1;

            return 0;
        }

        ql->append(val);

        Py_DECREF(itm);
    }

    Py_DECREF(iter);
 
    *sipCppPtr = ql;
 
    return sipGetState(sipTransferObj);
%End
};


%MappedType QList<qreal>
        /TypeHintIn="Iterable[float]", TypeHintOut="List[float]",
        TypeHintValue="[]"/
{
%TypeHeaderCode
#include <qlist.h>
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

    QList<qreal> *ql = new QList<qreal>;
 
    for (Py_ssize_t i = 0; ; ++i)
    {
        PyErr_Clear();
        PyObject *itm = PyIter_Next(iter);

        if (!itm)
        {
            if (PyErr_Occurred())
            {
                delete ql;
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
                    "index %zd has type '%s' but 'float' is expected", i,
                    sipPyTypeName(Py_TYPE(itm)));

            Py_DECREF(itm);
            delete ql;
            Py_DECREF(iter);
            *sipIsErr = 1;

            return 0;
        }

        ql->append(val);

        Py_DECREF(itm);
    }

    Py_DECREF(iter);
 
    *sipCppPtr = ql;
 
    return sipGetState(sipTransferObj);
%End
};


%MappedType QList<Qt::DayOfWeek>
        /TypeHintIn="Iterable[Qt.DayOfWeek]", TypeHintOut="List[Qt.DayOfWeek]",
        TypeHintValue="[]"/
{
%TypeHeaderCode
#include <Qt>
%End

%ConvertFromTypeCode
    PyObject *l = PyList_New(sipCpp->size());

    if (!l)
        return 0;

    for (int i = 0; i < sipCpp->size(); ++i)
    {
        PyObject *eobj = sipConvertFromEnum(sipCpp->at(i),
                sipType_Qt_DayOfWeek);

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

    QList<Qt::DayOfWeek> *ql = new QList<Qt::DayOfWeek>;

    for (Py_ssize_t i = 0; ; ++i)
    {
        PyErr_Clear();
        PyObject *itm = PyIter_Next(iter);

        if (!itm)
        {
            if (PyErr_Occurred())
            {
                delete ql;
                Py_DECREF(iter);
                *sipIsErr = 1;

                return 0;
            }

            break;
        }

        int v = sipConvertToEnum(itm, sipType_Qt_DayOfWeek);

        if (PyErr_Occurred())
        {
            PyErr_Format(PyExc_TypeError,
                    "index %zd has type '%s' but 'Qt.DayOfWeek' is expected",
                    i, sipPyTypeName(Py_TYPE(itm)));

            Py_DECREF(itm);
            delete ql;
            Py_DECREF(iter);
            *sipIsErr = 1;

            return 0;
        }

        ql->append(static_cast<Qt::DayOfWeek>(v));

        Py_DECREF(itm);
    }

    Py_DECREF(iter);

    *sipCppPtr = ql;

    return sipGetState(sipTransferObj);
%End
};
