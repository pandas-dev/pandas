// This is the SIP interface definition for the majority of the QPair based
// mapped types.
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


%If (Qt_5_10_0 -)

template<_TYPE1_, _TYPE2_>
%MappedType QPair<_TYPE1_, _TYPE2_> /TypeHint="Tuple[_TYPE1_, _TYPE2_]"/
{
%TypeHeaderCode
#include <qpair.h>
%End

%ConvertFromTypeCode
    _TYPE1_ *first = new _TYPE1_(sipCpp->first);
    _TYPE2_ *second = new _TYPE2_(sipCpp->second);
    PyObject *t = sipBuildResult(NULL, "(NN)", first, sipType__TYPE1_,
            sipTransferObj, second, sipType__TYPE2_, sipTransferObj);

    if (!t)
    {
        delete first;
        delete second;

        return 0;
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

    if (len != 2)
    {
        // A negative length should only be an internal error so let the
        // original exception stand.
        if (len >= 0)
            PyErr_Format(PyExc_TypeError,
                    "sequence has %zd elements but 2 elements are expected",
                    len);

        *sipIsErr = 1;

        return 0;
    }

    PyObject *firstobj = PySequence_GetItem(sipPy, 0);

    if (!firstobj)
    {
        *sipIsErr = 1;

        return 0;
    }

    int firststate;
    _TYPE1_ *first = reinterpret_cast<_TYPE1_ *>(
            sipForceConvertToType(firstobj, sipType__TYPE1_, sipTransferObj,
                    SIP_NOT_NONE, &firststate, sipIsErr));
 
    if (*sipIsErr)
    {
        PyErr_Format(PyExc_TypeError,
                "the first element has type '%s' but '_TYPE1_' is expected",
                sipPyTypeName(Py_TYPE(firstobj)));

        return 0;
    }

    PyObject *secondobj = PySequence_GetItem(sipPy, 1);

    if (!secondobj)
    {
        sipReleaseType(first, sipType__TYPE1_, firststate);
        Py_DECREF(firstobj);
        *sipIsErr = 1;

        return 0;
    }

    int secondstate;
    _TYPE2_ *second = reinterpret_cast<_TYPE2_ *>(
            sipForceConvertToType(secondobj, sipType__TYPE2_, sipTransferObj,
                    SIP_NOT_NONE, &secondstate, sipIsErr));
 
    if (*sipIsErr)
    {
        PyErr_Format(PyExc_TypeError,
                "the second element has type '%s' but '_TYPE2_' is expected",
                sipPyTypeName(Py_TYPE(secondobj)));

        Py_DECREF(secondobj);
        sipReleaseType(first, sipType__TYPE1_, firststate);
        Py_DECREF(firstobj);
        *sipIsErr = 1;

        return 0;
    }

    *sipCppPtr = new QPair<_TYPE1_, _TYPE2_>(*first, *second);

    sipReleaseType(second, sipType__TYPE2_, secondstate);
    Py_DECREF(secondobj);
    sipReleaseType(first, sipType__TYPE1_, firststate);
    Py_DECREF(firstobj);
 
    return sipGetState(sipTransferObj);
%End
};

%End


template<_TYPE_>
%MappedType QPair<_TYPE_, int> /TypeHint="Tuple[_TYPE_, int]"/
{
%TypeHeaderCode
#include <qpair.h>
%End

%ConvertFromTypeCode
    _TYPE_ *first = new _TYPE_(sipCpp->first);
    PyObject *t = sipBuildResult(NULL, "(Ni)", first, sipType__TYPE_,
            sipTransferObj, sipCpp->second);

    if (!t)
    {
        delete first;

        return 0;
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

    if (len != 2)
    {
        // A negative length should only be an internal error so let the
        // original exception stand.
        if (len >= 0)
            PyErr_Format(PyExc_TypeError,
                    "sequence has %zd elements but 2 elements are expected",
                    len);

        *sipIsErr = 1;

        return 0;
    }

    PyObject *firstobj = PySequence_GetItem(sipPy, 0);

    if (!firstobj)
    {
        *sipIsErr = 1;

        return 0;
    }

    int firststate;
    _TYPE_ *first = reinterpret_cast<_TYPE_ *>(
            sipForceConvertToType(firstobj, sipType__TYPE_, sipTransferObj,
                    SIP_NOT_NONE, &firststate, sipIsErr));
 
    if (*sipIsErr)
    {
        PyErr_Format(PyExc_TypeError,
                "the first element has type '%s' but '_TYPE_' is expected",
                sipPyTypeName(Py_TYPE(firstobj)));

        return 0;
    }

    PyObject *secondobj = PySequence_GetItem(sipPy, 1);

    if (!secondobj)
    {
        sipReleaseType(first, sipType__TYPE_, firststate);
        Py_DECREF(firstobj);
        *sipIsErr = 1;

        return 0;
    }

    int second = sipLong_AsInt(secondobj);

    if (PyErr_Occurred())
    {
        if (PyErr_ExceptionMatches(PyExc_TypeError))
            PyErr_Format(PyExc_TypeError,
                    "the second element has type '%s' but 'int' is expected",
                    sipPyTypeName(Py_TYPE(secondobj)));

        Py_DECREF(secondobj);
        sipReleaseType(first, sipType__TYPE_, firststate);
        Py_DECREF(firstobj);
        *sipIsErr = 1;

        return 0;
    }

    *sipCppPtr = new QPair<_TYPE_, int>(*first, second);

    Py_DECREF(secondobj);
    sipReleaseType(first, sipType__TYPE_, firststate);
    Py_DECREF(firstobj);
 
    return sipGetState(sipTransferObj);
%End
};


%MappedType QPair<int, int> /TypeHint="Tuple[int, int]"/
{
%TypeHeaderCode
#include <qpair.h>
%End

%ConvertFromTypeCode
    return Py_BuildValue("(ii)", sipCpp->first, sipCpp->second);
%End

%ConvertToTypeCode
    if (!sipIsErr)
        return (PySequence_Check(sipPy)
#if PY_MAJOR_VERSION < 3
                && !PyString_Check(sipPy)
#endif
                && !PyUnicode_Check(sipPy));

    Py_ssize_t len = PySequence_Size(sipPy);

    if (len != 2)
    {
        // A negative length should only be an internal error so let the
        // original exception stand.
        if (len >= 0)
            PyErr_Format(PyExc_TypeError,
                    "sequence has %zd elements but 2 elements are expected",
                    len);

        *sipIsErr = 1;

        return 0;
    }

    PyObject *firstobj = PySequence_GetItem(sipPy, 0);

    if (!firstobj)
    {
        *sipIsErr = 1;

        return 0;
    }

    int first = sipLong_AsInt(firstobj);

    if (PyErr_Occurred())
    {
        if (PyErr_ExceptionMatches(PyExc_TypeError))
            PyErr_Format(PyExc_TypeError,
                    "the first element has type '%s' but 'int' is expected",
                    sipPyTypeName(Py_TYPE(firstobj)));

        *sipIsErr = 1;

        return 0;
    }

    PyObject *secondobj = PySequence_GetItem(sipPy, 1);

    if (!secondobj)
    {
        Py_DECREF(firstobj);
        *sipIsErr = 1;

        return 0;
    }

    int second = sipLong_AsInt(secondobj);

    if (PyErr_Occurred())
    {
        if (PyErr_ExceptionMatches(PyExc_TypeError))
            PyErr_Format(PyExc_TypeError,
                    "the second element has type '%s' but 'int' is expected",
                    sipPyTypeName(Py_TYPE(secondobj)));

        Py_DECREF(secondobj);
        Py_DECREF(firstobj);
        *sipIsErr = 1;

        return 0;
    }

    *sipCppPtr = new QPair<int, int>(first, second);

    Py_DECREF(secondobj);
    Py_DECREF(firstobj);
 
    return sipGetState(sipTransferObj);
%End
};
