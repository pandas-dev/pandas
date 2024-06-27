// This is the SIP interface definition for the majority of the QHash based
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


template<_TYPE1_, _TYPE2_>
%MappedType QHash<_TYPE1_, _TYPE2_>
        /TypeHint="Dict[_TYPE1_, _TYPE2_]", TypeHintValue="{}"/
{
%TypeHeaderCode
#include <qhash.h>
%End

%ConvertFromTypeCode
    PyObject *d = PyDict_New();

    if (!d)
        return 0;

    QHash<_TYPE1_, _TYPE2_>::const_iterator it = sipCpp->constBegin();
    QHash<_TYPE1_, _TYPE2_>::const_iterator end = sipCpp->constEnd();

    while (it != end)
    {
        _TYPE1_ *k = new _TYPE1_(it.key());
        PyObject *kobj = sipConvertFromNewType(k, sipType__TYPE1_,
                sipTransferObj);

        if (!kobj)
        {
            delete k;
            Py_DECREF(d);

            return 0;
        }

        _TYPE2_ *v = new _TYPE2_(it.value());
        PyObject *vobj = sipConvertFromNewType(v, sipType__TYPE2_,
                sipTransferObj);

        if (!vobj)
        {
            delete v;
            Py_DECREF(kobj);
            Py_DECREF(d);

            return 0;
        }

        int rc = PyDict_SetItem(d, kobj, vobj);

        Py_DECREF(vobj);
        Py_DECREF(kobj);

        if (rc < 0)
        {
            Py_DECREF(d);

            return 0;
        }

        ++it;
    }

    return d;
%End

%ConvertToTypeCode
    if (!sipIsErr)
        return PyDict_Check(sipPy);

    QHash<_TYPE1_, _TYPE2_> *qh = new QHash<_TYPE1_, _TYPE2_>;

    Py_ssize_t pos = 0;
    PyObject *kobj, *vobj;
 
    while (PyDict_Next(sipPy, &pos, &kobj, &vobj))
    {
        int kstate;
        _TYPE1_ *k = reinterpret_cast<_TYPE1_ *>(
                sipForceConvertToType(kobj, sipType__TYPE1_, sipTransferObj,
                        SIP_NOT_NONE, &kstate, sipIsErr));

        if (*sipIsErr)
        {
            PyErr_Format(PyExc_TypeError,
                    "a dict key has type '%s' but '_TYPE1_' is expected",
                    sipPyTypeName(Py_TYPE(kobj)));

            delete qh;

            return 0;
        }

        int vstate;
        _TYPE2_ *v = reinterpret_cast<_TYPE2_ *>(
                sipForceConvertToType(vobj, sipType__TYPE2_, sipTransferObj,
                        SIP_NOT_NONE, &vstate, sipIsErr));

        if (*sipIsErr)
        {
            PyErr_Format(PyExc_TypeError,
                    "a dict value has type '%s' but '_TYPE2_' is expected",
                    sipPyTypeName(Py_TYPE(vobj)));

            sipReleaseType(k, sipType__TYPE1_, kstate);
            delete qh;

            return 0;
        }

        qh->insert(*k, *v);

        sipReleaseType(v, sipType__TYPE2_, vstate);
        sipReleaseType(k, sipType__TYPE1_, kstate);
    }
 
    *sipCppPtr = qh;
 
    return sipGetState(sipTransferObj);
%End
};


%If (Qt_5_4_0 -)

// This is only needed for QtWebChannel but is sufficiently generic that we
// include it here.

template<_TYPE1_, _TYPE2_ *>
%MappedType QHash<_TYPE1_, _TYPE2_ *>
        /TypeHint="Dict[_TYPE1_, _TYPE2_]", TypeHintValue="{}"/
{
%TypeHeaderCode
#include <qhash.h>
%End

%ConvertFromTypeCode
    int gc_enabled = sipEnableGC(0);
    PyObject *d = PyDict_New();

    if (d)
    {
        QHash<_TYPE1_, _TYPE2_ *>::const_iterator it = sipCpp->constBegin();
        QHash<_TYPE1_, _TYPE2_ *>::const_iterator end = sipCpp->constEnd();

        while (it != end)
        {
            _TYPE1_ *k = new _TYPE1_(it.key());
            PyObject *kobj = sipConvertFromNewType(k, sipType__TYPE1_,
                    sipTransferObj);

            if (!kobj)
            {
                delete k;
                Py_DECREF(d);
                d = 0;

                break;
            }

            _TYPE2_ *v = it.value();
            PyObject *vobj = sipConvertFromNewType(v, sipType__TYPE2_,
                    sipTransferObj);

            if (!vobj)
            {
                Py_DECREF(kobj);
                Py_DECREF(d);
                d = 0;

                break;
            }

            int rc = PyDict_SetItem(d, kobj, vobj);

            Py_DECREF(vobj);
            Py_DECREF(kobj);

            if (rc < 0)
            {
                Py_DECREF(d);
                d = 0;

                break;
            }

            ++it;
        }
    }

    sipEnableGC(gc_enabled);

    return d;
%End

%ConvertToTypeCode
    if (!sipIsErr)
        return PyDict_Check(sipPy);

    QHash<_TYPE1_, _TYPE2_ *> *qh = new QHash<_TYPE1_, _TYPE2_ *>;

    Py_ssize_t pos = 0;
    PyObject *kobj, *vobj;
 
    while (PyDict_Next(sipPy, &pos, &kobj, &vobj))
    {
        int kstate;
        _TYPE1_ *k = reinterpret_cast<_TYPE1_ *>(
                sipForceConvertToType(kobj, sipType__TYPE1_, sipTransferObj,
                        SIP_NOT_NONE, &kstate, sipIsErr));

        if (*sipIsErr)
        {
            PyErr_Format(PyExc_TypeError,
                    "a dict key has type '%s' but '_TYPE1_' is expected",
                    sipPyTypeName(Py_TYPE(kobj)));

            delete qh;

            return 0;
        }

        _TYPE2_ *v = reinterpret_cast<_TYPE2_ *>(
                sipForceConvertToType(vobj, sipType__TYPE2_, sipTransferObj, 0,
                        0, sipIsErr));

        if (*sipIsErr)
        {
            PyErr_Format(PyExc_TypeError,
                    "a dict value has type '%s' but '_TYPE2_' is expected",
                    sipPyTypeName(Py_TYPE(vobj)));

            sipReleaseType(k, sipType__TYPE1_, kstate);
            delete qh;

            return 0;
        }

        qh->insert(*k, v);

        sipReleaseType(k, sipType__TYPE1_, kstate);
    }
 
    *sipCppPtr = qh;
 
    return sipGetState(sipTransferObj);
%End
};

%End


template<int, _TYPE_>
%MappedType QHash<int, _TYPE_>
        /TypeHint="Dict[int, _TYPE_]", TypeHintValue="{}"/
{
%TypeHeaderCode
#include <qhash.h>
%End

%ConvertFromTypeCode
    PyObject *d = PyDict_New();

    if (!d)
        return 0;

    QHash<int, _TYPE_>::const_iterator it = sipCpp->constBegin();
    QHash<int, _TYPE_>::const_iterator end = sipCpp->constEnd();

    while (it != end)
    {
        PyObject *kobj = SIPLong_FromLong(it.key());

        if (!kobj)
        {
            Py_DECREF(d);

            return 0;
        }

        _TYPE_ *v = new _TYPE_(it.value());
        PyObject *vobj = sipConvertFromNewType(v, sipType__TYPE_,
                sipTransferObj);

        if (!vobj)
        {
            delete v;
            Py_DECREF(kobj);
            Py_DECREF(d);

            return 0;
        }

        int rc = PyDict_SetItem(d, kobj, vobj);

        Py_DECREF(vobj);
        Py_DECREF(kobj);

        if (rc < 0)
        {
            Py_DECREF(d);

            return 0;
        }

        ++it;
    }

    return d;
%End

%ConvertToTypeCode
    if (!sipIsErr)
        return PyDict_Check(sipPy);

    QHash<int, _TYPE_> *qh = new QHash<int, _TYPE_>;

    Py_ssize_t pos = 0;
    PyObject *kobj, *vobj;
 
    while (PyDict_Next(sipPy, &pos, &kobj, &vobj))
    {
        int k = sipLong_AsInt(kobj);

        if (PyErr_Occurred())
        {
            if (PyErr_ExceptionMatches(PyExc_TypeError))
                PyErr_Format(PyExc_TypeError,
                        "a dict key has type '%s' but 'int' is expected",
                        sipPyTypeName(Py_TYPE(kobj)));

            delete qh;
            *sipIsErr = 1;

            return 0;
        }

        int vstate;
        _TYPE_ *v = reinterpret_cast<_TYPE_ *>(
                sipForceConvertToType(vobj, sipType__TYPE_, sipTransferObj,
                        SIP_NOT_NONE, &vstate, sipIsErr));

        if (*sipIsErr)
        {
            PyErr_Format(PyExc_TypeError,
                    "a dict value has type '%s' but '_TYPE_' is expected",
                    sipPyTypeName(Py_TYPE(vobj)));

            delete qh;

            return 0;
        }

        qh->insert(k, *v);

        sipReleaseType(v, sipType__TYPE_, vstate);
    }
 
    *sipCppPtr = qh;
 
    return sipGetState(sipTransferObj);
%End
};


%If (Qt_5_12_0 -)

template<quint16, _TYPE_>
%MappedType QHash<quint16, _TYPE_>
        /TypeHint="Dict[int, _TYPE_]", TypeHintValue="{}"/
{
%TypeHeaderCode
#include <qhash.h>
%End

%ConvertFromTypeCode
    PyObject *d = PyDict_New();

    if (!d)
        return 0;

    QHash<quint16, _TYPE_>::const_iterator it = sipCpp->constBegin();
    QHash<quint16, _TYPE_>::const_iterator end = sipCpp->constEnd();

    while (it != end)
    {
        PyObject *kobj = SIPLong_FromLong(it.key());

        if (!kobj)
        {
            Py_DECREF(d);

            return 0;
        }

        _TYPE_ *v = new _TYPE_(it.value());
        PyObject *vobj = sipConvertFromNewType(v, sipType__TYPE_,
                sipTransferObj);

        if (!vobj)
        {
            delete v;
            Py_DECREF(kobj);
            Py_DECREF(d);

            return 0;
        }

        int rc = PyDict_SetItem(d, kobj, vobj);

        Py_DECREF(vobj);
        Py_DECREF(kobj);

        if (rc < 0)
        {
            Py_DECREF(d);

            return 0;
        }

        ++it;
    }

    return d;
%End

%ConvertToTypeCode
    if (!sipIsErr)
        return PyDict_Check(sipPy);

    QHash<quint16, _TYPE_> *qh = new QHash<quint16, _TYPE_>;

    Py_ssize_t pos = 0;
    PyObject *kobj, *vobj;
 
    while (PyDict_Next(sipPy, &pos, &kobj, &vobj))
    {
        quint16 k = sipLong_AsUnsignedShort(kobj);

        if (PyErr_Occurred())
        {
            if (PyErr_ExceptionMatches(PyExc_TypeError))
                PyErr_Format(PyExc_TypeError,
                        "a dict key has type '%s' but 'int' is expected",
                        sipPyTypeName(Py_TYPE(kobj)));

            delete qh;
            *sipIsErr = 1;

            return 0;
        }

        int vstate;
        _TYPE_ *v = reinterpret_cast<_TYPE_ *>(
                sipForceConvertToType(vobj, sipType__TYPE_, sipTransferObj,
                        SIP_NOT_NONE, &vstate, sipIsErr));

        if (*sipIsErr)
        {
            PyErr_Format(PyExc_TypeError,
                    "a dict value has type '%s' but '_TYPE_' is expected",
                    sipPyTypeName(Py_TYPE(vobj)));

            delete qh;

            return 0;
        }

        qh->insert(k, *v);

        sipReleaseType(v, sipType__TYPE_, vstate);
    }
 
    *sipCppPtr = qh;
 
    return sipGetState(sipTransferObj);
%End
};

%End
