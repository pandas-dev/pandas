// This is the SIP interface definition for the QJsonObject mapped type.
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


%MappedType QJsonObject
        /TypeHint="Dict[QString, QJsonValue]", TypeHintValue="{}"/
{
%TypeHeaderCode
#include <qjsonobject.h>
%End

%ConvertFromTypeCode
    PyObject *d = PyDict_New();

    if (!d)
        return 0;

    QJsonObject::const_iterator it = sipCpp->constBegin();
    QJsonObject::const_iterator end = sipCpp->constEnd();

    while (it != end)
    {
        QString *k = new QString(it.key());
        PyObject *kobj = sipConvertFromNewType(k, sipType_QString,
                sipTransferObj);

        if (!kobj)
        {
            delete k;
            Py_DECREF(d);

            return 0;
        }

        QJsonValue *v = new QJsonValue(it.value());
        PyObject *vobj = sipConvertFromNewType(v, sipType_QJsonValue,
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

    QJsonObject *jo = new QJsonObject;

    Py_ssize_t pos = 0;
    PyObject *kobj, *vobj;
 
    while (PyDict_Next(sipPy, &pos, &kobj, &vobj))
    {
        int kstate;
        QString *k = reinterpret_cast<QString *>(
                sipForceConvertToType(kobj, sipType_QString, sipTransferObj,
                        SIP_NOT_NONE, &kstate, sipIsErr));

        if (*sipIsErr)
        {
            PyErr_Format(PyExc_TypeError,
                    "a key has type '%s' but 'str' is expected",
                    sipPyTypeName(Py_TYPE(kobj)));

            delete jo;

            return 0;
        }

        int vstate;
        QJsonValue *v = reinterpret_cast<QJsonValue *>(
                sipForceConvertToType(vobj, sipType_QJsonValue, sipTransferObj,
                        SIP_NOT_NONE, &vstate, sipIsErr));

        if (*sipIsErr)
        {
            PyErr_Format(PyExc_TypeError,
                    "a value has type '%s' but 'QJsonValue' is expected",
                    sipPyTypeName(Py_TYPE(vobj)));

            sipReleaseType(k, sipType_QString, kstate);
            delete jo;

            return 0;
        }

        jo->insert(*k, *v);

        sipReleaseType(v, sipType_QJsonValue, vstate);
        sipReleaseType(k, sipType_QString, kstate);
    }
 
    *sipCppPtr = jo;
 
    return sipGetState(sipTransferObj);
%End
};
