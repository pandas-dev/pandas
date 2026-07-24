// This is the SIP interface definition for QWinJumpListItem.
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

class QWinJumpListItem /Supertype=sip.wrapper/
{
%TypeHeaderCode
#include <qwinjumplistitem.h>
%End

public:
    enum Type
    {
        Destination,
        Link,
        Separator
    };

    explicit QWinJumpListItem(Type type);
    ~QWinJumpListItem();

    void setType(Type type);
    Type type() const;
    void setFilePath(const QString &filePath);
    QString filePath() const;
    void setWorkingDirectory(const QString &workingDirectory);
    QString workingDirectory() const;
    void setIcon(const QIcon &icon);
    QIcon icon() const;
    void setTitle(const QString &title);
    QString title() const;
    void setDescription(const QString &description);
    QString description() const;
    void setArguments(const QStringList &arguments);
    QStringList arguments() const;

private:
    QWinJumpListItem(const QWinJumpListItem &);
};

%End
