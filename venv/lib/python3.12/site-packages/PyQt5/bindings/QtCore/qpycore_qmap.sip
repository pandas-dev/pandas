// This is the SIP interface definition for the majority of the QMap and
// QMultiMap based mapped types.
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


template<_TYPE1_, _TYPE2_>
%MappedType QMap<_TYPE1_, _TYPE2_>
        /TypeHint="Dict[_TYPE1_, _TYPE2_]", TypeHintValue="{}"/
{
%TypeHeaderCode
#include <qmap.h>
%End

%ConvertFromTypeCode
    PyObject *d = PyDict_New();

    if (!d)
        return 0;

    QMap<_TYPE1_, _TYPE2_>::const_iterator it = sipCpp->constBegin();
    QMap<_TYPE1_, _TYPE2_>::const_iterator end = sipCpp->constEnd();

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

    QMap<_TYPE1_, _TYPE2_> *qm = new QMap<_TYPE1_, _TYPE2_>;

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

            delete qm;

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
            delete qm;

            return 0;
        }

        qm->insert(*k, *v);

        sipReleaseType(v, sipType__TYPE2_, vstate);
        sipReleaseType(k, sipType__TYPE1_, kstate);
    }
 
    *sipCppPtr = qm;
 
    return sipGetState(sipTransferObj);
%End
};


template<int, _TYPE_>
%MappedType QMap<int, _TYPE_>
        /TypeHint="Dict[int, _TYPE_]", TypeHintValue="{}"/
{
%TypeHeaderCode
#include <qmap.h>
%End

%ConvertFromTypeCode
    PyObject *d = PyDict_New();

    if (!d)
        return 0;

    QMap<int, _TYPE_>::const_iterator it = sipCpp->constBegin();
    QMap<int, _TYPE_>::const_iterator end = sipCpp->constEnd();

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

    QMap<int, _TYPE_> *qm = new QMap<int, _TYPE_>;

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

            delete qm;
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

            delete qm;

            return 0;
        }

        qm->insert(k, *v);

        sipReleaseType(v, sipType__TYPE_, vstate);
    }
 
    *sipCppPtr = qm;
 
    return sipGetState(sipTransferObj);
%End
};


template<_TYPE1_, _TYPE2_>
%MappedType QMultiMap<_TYPE1_, _TYPE2_>
        /TypeHintOut="Dict[_TYPE1_, List[_TYPE2_]]", TypeHintValue="{}"/
{
%TypeHeaderCode
#include <qmap.h>
%End

%ConvertFromTypeCode
    PyObject *d = PyDict_New();

    if (!d)
        return 0;

    QList<_TYPE1_> keys = sipCpp->keys();
    QList<_TYPE1_>::const_iterator kit = keys.constBegin();
    QList<_TYPE1_>::const_iterator kit_end = keys.constEnd();

    while (kit != kit_end)
    {
        _TYPE1_ *k = new _TYPE1_(*kit);
        PyObject *kobj = sipConvertFromNewType(k, sipType__TYPE1_,
                sipTransferObj);

        if (!kobj)
        {
            delete k;
            Py_DECREF(d);

            return 0;
        }

        // Create a Python list as the dictionary value.
        QList<_TYPE2_> values = sipCpp->values(*kit);
        PyObject *vobj = PyList_New(values.count());

        if (!vobj)
        {
            Py_DECREF(kobj);
            Py_DECREF(d);

            return 0;
        }

        QList<_TYPE2_>::const_iterator vit = values.constBegin();
        QList<_TYPE2_>::const_iterator vit_end = values.constEnd();

        for (int i = 0; vit != vit_end; ++i)
        {
            _TYPE2_ *sv = new _TYPE2_(*vit);
            PyObject *svobj = sipConvertFromNewType(sv, sipType__TYPE2_,
                    sipTransferObj);

            if (!svobj)
            {
                delete sv;
                Py_DECREF(vobj);
                Py_DECREF(kobj);
                Py_DECREF(d);

                return 0;
            }

            PyList_SetItem(vobj, i, svobj);

            ++vit;
        }

        int rc = PyDict_SetItem(d, kobj, vobj);

        Py_DECREF(vobj);
        Py_DECREF(kobj);

        if (rc < 0)
        {
            Py_DECREF(d);

            return 0;
        }

        ++kit;
    }

    return d;
%End

%ConvertToTypeCode
    if (!sipIsErr)
        return PyDict_Check(sipPy);

    // Note that PyQt v5.1 contains an unused implementation that can be
    // restored if needed (although it will need updating to accept an iterable
    // rather than just a list of values).
    PyErr_SetString(PyExc_NotImplementedError,
            "converting to QMultiMap<_TYPE1_, _TYPE2_> is unsupported");

    *sipIsErr = 1;

    return 0;
%End
};
