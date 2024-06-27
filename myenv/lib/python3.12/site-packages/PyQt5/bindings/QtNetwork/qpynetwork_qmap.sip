// This is the SIP interface definition for the QMap and QMultiMap based mapped
// types specific to the QtNetwork module.
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


%If (PyQt_SSL)

%MappedType QMultiMap<QSsl::AlternativeNameEntryType, QString>
        /TypeHintOut="Dict[QSsl.AlternativeNameEntryType, List[QString]]",
        TypeHintValue="[]"/
{
%TypeHeaderCode
#include <qmap.h>
#include <qssl.h>
%End

%ConvertFromTypeCode
    // Get the enum member objects that are the dictionary keys.
    static PyObject *email_entry = 0;

    if (!email_entry)
    {
        email_entry = PyObject_GetAttrString(
                (PyObject *)sipTypeAsPyTypeObject(sipType_QSsl), "EmailEntry");

        if (!email_entry)
            return 0;
    }

    static PyObject *dns_entry = 0;

    if (!dns_entry)
    {
        dns_entry = PyObject_GetAttrString(
                (PyObject *)sipTypeAsPyTypeObject(sipType_QSsl), "DnsEntry");

        if (!dns_entry)
            return 0;
    }

#if QT_VERSION >= 0x050d00
    static PyObject *ip_address_entry = 0;

    if (!ip_address_entry)
    {
        ip_address_entry = PyObject_GetAttrString(
                (PyObject *)sipTypeAsPyTypeObject(sipType_QSsl), "IpAddressEntry");

        if (!ip_address_entry)
            return 0;
    }
#endif

    // Create the dictionary.
    PyObject *d = PyDict_New();

    if (!d)
        return 0;

    QList<QString> vl;

    // Handle the Qssl::EmailEntry key.
    vl = sipCpp->values(QSsl::EmailEntry);

    if (!vl.isEmpty())
    {
        PyObject *vlobj = PyList_New(vl.count());

        if (!vlobj)
        {
            Py_DECREF(d);
            return 0;
        }

        int rc = PyDict_SetItem(d, email_entry, vlobj);

        Py_DECREF(vlobj);

        if (rc < 0)
        {
            Py_DECREF(d);
            return 0;
        }

        for (int i = 0; i < vl.count(); ++i)
        {
            QString *s = new QString(vl.at(i));
            PyObject *vobj = sipConvertFromNewType(s, sipType_QString,
                    sipTransferObj);

            if (!vobj)
            {
                delete s;
                Py_DECREF(d);
                return 0;
            }

            PyList_SetItem(vlobj, i, vobj);
        }
    }

    // Handle the Qssl::DnsEntry key.
    vl = sipCpp->values(QSsl::DnsEntry);

    if (!vl.isEmpty())
    {
        PyObject *vlobj = PyList_New(vl.count());

        if (!vlobj)
        {
            Py_DECREF(d);
            return 0;
        }

        int rc = PyDict_SetItem(d, dns_entry, vlobj);

        Py_DECREF(vlobj);

        if (rc < 0)
        {
            Py_DECREF(d);
            return 0;
        }

        for (int i = 0; i < vl.count(); ++i)
        {
            QString *s = new QString(vl.at(i));
            PyObject *vobj = sipConvertFromNewType(s, sipType_QString,
                    sipTransferObj);

            if (!vobj)
            {
                delete s;
                Py_DECREF(d);
                return 0;
            }

            PyList_SetItem(vlobj, i, vobj);
        }
    }

#if QT_VERSION >= 0x050d00
    // Handle the Qssl::IpAddressEntry key.
    vl = sipCpp->values(QSsl::IpAddressEntry);

    if (!vl.isEmpty())
    {
        PyObject *vlobj = PyList_New(vl.count());

        if (!vlobj)
        {
            Py_DECREF(d);
            return 0;
        }

        int rc = PyDict_SetItem(d, ip_address_entry, vlobj);

        Py_DECREF(vlobj);

        if (rc < 0)
        {
            Py_DECREF(d);
            return 0;
        }

        for (int i = 0; i < vl.count(); ++i)
        {
            QString *s = new QString(vl.at(i));
            PyObject *vobj = sipConvertFromNewType(s, sipType_QString,
                    sipTransferObj);

            if (!vobj)
            {
                delete s;
                Py_DECREF(d);
                return 0;
            }

            PyList_SetItem(vlobj, i, vobj);
        }
    }
#endif

    return d;
%End

%ConvertToTypeCode
    if (!sipIsErr)
        return PyDict_Check(sipPy);

    PyErr_SetString(PyExc_NotImplementedError,
            "converting to QMultiMap<QSsl::AlternativeNameEntryType, QString> is unsupported");

    return 0;
%End
};

%End
