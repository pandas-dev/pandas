// This is the SIP specification of the QPyDBusReply class.
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


class QPyDBusReply /PyName=QDBusReply/
{
%TypeHeaderCode
#include <qpydbusreply.h>
%End

public:
    QPyDBusReply(const QDBusMessage &reply) /HoldGIL/;
    QPyDBusReply(const QDBusPendingCall &call) /HoldGIL/;
    QPyDBusReply(const QDBusError &error);
    QPyDBusReply(const QPyDBusReply &other) /HoldGIL/;
    ~QPyDBusReply() /HoldGIL/;

    const QDBusError &error() const /HoldGIL/;
    bool isValid() const /HoldGIL/;
    SIP_PYOBJECT value(SIP_PYOBJECT type /TypeHintValue="None"/ = 0) const /HoldGIL/;
};


template<TYPE>
%MappedType QDBusReply<TYPE> /TypeHint="QDBusReply"/
{
%TypeHeaderCode
#include <qdbusreply.h>
#include <qpydbusreply.h>
%End

%ConvertFromTypeCode
    PyObject *value_obj;

    if (sipCpp->isValid())
    {
        // Convert the value to a Python object.
        TYPE *value = new TYPE(sipCpp->value());

        if ((value_obj = sipConvertFromNewType(value, sipType_TYPE, NULL)) == NULL)
        {
            delete value;
            return NULL;
        }
    }
    else
    {
        value_obj = 0;
    }

    QPyDBusReply *reply = new QPyDBusReply(value_obj,
            sipCpp->isValid(), sipCpp->error());

    PyObject *reply_obj = sipConvertFromNewType(reply, sipType_QPyDBusReply, sipTransferObj);

    if (reply_obj == NULL)
    {
        delete reply;
        return NULL;
    }

    return reply_obj;
%End

%ConvertToTypeCode
    // We never create a QDBusReply from Python.
    return 0;
%End
};


%MappedType QDBusReply<void> /TypeHint="QDBusReply"/
{
%TypeHeaderCode
#include <qdbusreply.h>
#include <qpydbusreply.h>
%End

%ConvertFromTypeCode
    Py_INCREF(Py_None);
    QPyDBusReply *reply = new QPyDBusReply(Py_None,
            sipCpp->isValid(), sipCpp->error());

    PyObject *reply_obj = sipConvertFromNewType(reply, sipType_QPyDBusReply, sipTransferObj);

    if (reply_obj == NULL)
    {
        delete reply;
        return NULL;
    }

    return reply_obj;
%End

%ConvertToTypeCode
    // We never create a QDBusReply from Python.
    return 0;
%End
};


%MappedType QDBusReply<bool> /TypeHint="QDBusReply"/
{
%TypeHeaderCode
#include <qdbusreply.h>
#include <qpydbusreply.h>
%End

%ConvertFromTypeCode
    PyObject *value_obj;

    if (sipCpp->isValid())
    {
        if ((value_obj = PyBool_FromLong(sipCpp->value())) == NULL)
            return NULL;
    }
    else
    {
        value_obj = 0;
    }

    QPyDBusReply *reply = new QPyDBusReply(value_obj,
            sipCpp->isValid(), sipCpp->error());

    PyObject *reply_obj = sipConvertFromNewType(reply, sipType_QPyDBusReply, sipTransferObj);

    if (reply_obj == NULL)
    {
        delete reply;
        return NULL;
    }

    return reply_obj;
%End

%ConvertToTypeCode
    // We never create a QDBusReply from Python.
    return 0;
%End
};


%MappedType QDBusReply<unsigned> /TypeHint="QDBusReply"/
{
%TypeHeaderCode
#include <qdbusreply.h>
#include <qpydbusreply.h>
%End

%ConvertFromTypeCode
    PyObject *value_obj;

    if (sipCpp->isValid())
    {
        if ((value_obj = PyLong_FromUnsignedLong(sipCpp->value())) == NULL)
            return NULL;
    }
    else
    {
        value_obj = 0;
    }

    QPyDBusReply *reply = new QPyDBusReply(value_obj,
            sipCpp->isValid(), sipCpp->error());

    PyObject *reply_obj = sipConvertFromNewType(reply, sipType_QPyDBusReply, sipTransferObj);

    if (reply_obj == NULL)
    {
        delete reply;
        return NULL;
    }

    return reply_obj;
%End

%ConvertToTypeCode
    // We never create a QDBusReply from Python.
    return 0;
%End
};


%MappedType QDBusReply<QDBusConnectionInterface::RegisterServiceReply> /TypeHint="QDBusReply"/
{
%TypeHeaderCode
#include <qdbusreply.h>
#include <qpydbusreply.h>
%End

%ConvertFromTypeCode
    PyObject *value_obj;

    if (sipCpp->isValid())
    {
        if ((value_obj = sipConvertFromEnum(sipCpp->value(), sipType_QDBusConnectionInterface_RegisterServiceReply)) == NULL)
            return NULL;
    }
    else
    {
        value_obj = 0;
    }

    QPyDBusReply *reply = new QPyDBusReply(value_obj,
            sipCpp->isValid(), sipCpp->error());

    PyObject *reply_obj = sipConvertFromNewType(reply, sipType_QPyDBusReply, sipTransferObj);

    if (reply_obj == NULL)
    {
        delete reply;
        return NULL;
    }

    return reply_obj;
%End

%ConvertToTypeCode
    // We never create a QDBusReply from Python.
    return 0;
%End
};
