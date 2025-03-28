// This is the SIP interface definition for QMacPasteboardMime.
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


%If (Qt_5_2_0 -)

%If (PyQt_MacOSXOnly)

class QMacPasteboardMime
{
%TypeHeaderCode
#include <qmacpasteboardmime.h>
%End

public:
    enum QMacPasteboardMimeType
    {
        MIME_DND,
        MIME_CLIP,
        MIME_QT_CONVERTOR,
        MIME_QT3_CONVERTOR,
        MIME_ALL,
    };

    explicit QMacPasteboardMime(char t /PyInt/);
    virtual ~QMacPasteboardMime();

    virtual QString convertorName() = 0;
    virtual bool canConvert(const QString &mime, QString flav) = 0;
    virtual QString mimeFor(QString flav) = 0;
    virtual QString flavorFor(const QString &mime) = 0;
    virtual QVariant convertToMime(const QString &mime, QList<QByteArray> data,
            QString flav) = 0;
    virtual QList<QByteArray> convertFromMime(const QString &mime,
            QVariant data, QString flav) = 0;
    virtual int count(QMimeData *mimeData);
};


%ModuleHeaderCode
#include <qmacpasteboardmime.h>
%End

void qRegisterDraggedTypes(const QStringList &types);

%End

%End
