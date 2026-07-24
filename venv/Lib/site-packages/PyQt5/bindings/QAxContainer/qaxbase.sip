// This is the SIP interface definition for QAxBase.
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


class QAxBase /Abstract, PyQtNoQMetaObject/
{
%TypeHeaderCode
#include <qaxbase.h>
%End

public:
    //QAxBase(IUnknown *iface = 0);
    virtual ~QAxBase();

    QString control() const;

    //long queryInterface(const QUuid &, void **) const;

    // Note that the order of these overloads is significant.
    QVariant dynamicCall(const char *, QList<QVariant> & /GetWrapper/);
%MethodCode
        Py_BEGIN_ALLOW_THREADS
        sipRes = new QVariant(sipCpp->dynamicCall(a0, *a1));
        Py_END_ALLOW_THREADS

        // Update the input list with the (possibly) new values.
        for (Py_ssize_t i = 0; i < PyList_Size(a1Wrapper); ++i)
        {
            QVariant *v = new QVariant(a1->at(i));
            PyObject *v_obj = sipConvertFromNewType(v, sipType_QVariant, NULL);

            if (!v_obj)
            {
                delete v;
                sipIsErr = 1;
                break;
            }

            if (PyList_SetItem(a1Wrapper, i, v_obj) < 0)
            {
                Py_DECREF(v_obj);
                sipIsErr = 1;
                break;
            }
        }
%End

    QVariant dynamicCall(const char *,
            const QVariant &value1 = QVariant(),
            const QVariant &value2 = QVariant(),
            const QVariant &value3 = QVariant(),
            const QVariant &value4 = QVariant(),
            const QVariant &value5 = QVariant(),
            const QVariant &value6 = QVariant(),
            const QVariant &value7 = QVariant(),
            const QVariant &value8 = QVariant());

    // Note that the order of these overloads is significant.
    QAxObject *querySubObject(const char *, QList<QVariant> & /GetWrapper/);
%MethodCode
        Py_BEGIN_ALLOW_THREADS
        sipRes = sipCpp->querySubObject(a0, *a1);
        Py_END_ALLOW_THREADS

        // Update the input list with the (possibly) new values.
        for (Py_ssize_t i = 0; i < PyList_Size(a1Wrapper); ++i)
        {
            QVariant *v = new QVariant(a1->at(i));
            PyObject *v_obj = sipConvertFromNewType(v, sipType_QVariant, NULL);

            if (!v_obj)
            {
                delete v;
                sipIsErr = 1;
                break;
            }

            if (PyList_SetItem(a1Wrapper, i, v_obj) < 0)
            {
                Py_DECREF(v_obj);
                sipIsErr = 1;
                break;
            }
        }
%End

    QAxObject *querySubObject(const char *,
            const QVariant &value1 = QVariant(),
            const QVariant &value2 = QVariant(),
            const QVariant &value3 = QVariant(),
            const QVariant &value4 = QVariant(),
            const QVariant &value5 = QVariant(),
            const QVariant &value6 = QVariant(),
            const QVariant &value7 = QVariant(),
            const QVariant &value8 = QVariant());

    // SIP has a bug triggered by a template definition being the subject of
    // multiple typedefs.  It only really matters when building everything as
    // one big module (the code that implements the type is duplicated in
    // other cases).  Until it is fixed we just avoid the problematic typedef.
    //typedef QMap<QString, QVariant> PropertyBag;
    //PropertyBag propertyBag() const;
    //void setPropertyBag(const PropertyBag &);
    QVariantMap propertyBag() const;
    void setPropertyBag(const QVariantMap &);

    QString generateDocumentation();

    virtual bool propertyWritable(const char *) const;
    virtual void setPropertyWritable(const char *, bool);

    bool isNull() const;

    QStringList verbs() const;

    QVariant asVariant() const;

signals:
    void signal(const QString &, int, void *);
    void propertyChanged(const QString &);
    void exception(int, const QString &, const QString &, const QString &);

public:
    virtual void clear();
    bool setControl(const QString &);

    void disableMetaObject();
    void disableClassInfo();
    void disableEventSink();

%If (Qt_5_13_0 -)
    unsigned long classContext() const;
    void setClassContext(unsigned long classContext);
%End

protected:
    //virtual bool initialize(IUnknown** ptr);
    //bool initializeRemote(IUnknown** ptr);
    //bool initializeLicensed(IUnknown** ptr);
    //bool initializeActive(IUnknown** ptr);
    //bool initializeFromFile(IUnknown** ptr);
};
