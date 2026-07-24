// This is the SIP interface definition for the QHash based mapped types
// specific to the QtNetwork module.
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


%MappedType QHash<QNetworkRequest::Attribute, QVariant>
        /TypeHint="Dict[QNetworkRequest.Attribute, QVariant]",
        TypeHintValue="{}"/
{
%TypeHeaderCode
#include <qhash.h>
#include <qnetworkrequest.h>
#include <qvariant.h>
%End

%ConvertFromTypeCode
    PyObject *d = PyDict_New();

    if (!d)
        return 0;

    QHash<QNetworkRequest::Attribute, QVariant>::const_iterator it = sipCpp->constBegin();
    QHash<QNetworkRequest::Attribute, QVariant>::const_iterator end = sipCpp->constEnd();

    while (it != end)
    {
        PyObject *kobj = sipConvertFromEnum(it.key(),
                sipType_QNetworkRequest_Attribute);

        if (!kobj)
        {
            Py_DECREF(d);

            return 0;
        }

        QVariant *v = new QVariant(it.value());
        PyObject *vobj = sipConvertFromNewType(v, sipType_QVariant,
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

    QHash<QNetworkRequest::Attribute, QVariant> *qh = new QHash<QNetworkRequest::Attribute, QVariant>;

    Py_ssize_t pos = 0;
    PyObject *kobj, *vobj;
 
    while (PyDict_Next(sipPy, &pos, &kobj, &vobj))
    {
        int k = sipConvertToEnum(kobj, sipType_QNetworkRequest_Attribute);

        if (PyErr_Occurred())
        {
            PyErr_Format(PyExc_TypeError,
                    "a key has type '%s' but 'QNetworkRequest.Attribute' is expected",
                    sipPyTypeName(Py_TYPE(kobj)));

            delete qh;
            *sipIsErr = 1;

            return 0;
        }

        int vstate;
        QVariant *v = reinterpret_cast<QVariant *>(
                sipForceConvertToType(vobj, sipType_QVariant, sipTransferObj,
                        SIP_NOT_NONE, &vstate, sipIsErr));

        if (*sipIsErr)
        {
            // Any error must be internal, so leave the exception as it is.

            delete qh;

            return 0;
        }

        qh->insert(static_cast<QNetworkRequest::Attribute>(k), *v);

        sipReleaseType(v, sipType_QVariant, vstate);
    }
 
    *sipCppPtr = qh;
 
    return sipGetState(sipTransferObj);
%End
};
