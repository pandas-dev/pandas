// This is the SIP interface definition for QWinJumpList.
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

class QWinJumpList : public QObject
{
%TypeHeaderCode
#include <qwinjumplist.h>
%End

%ConvertToSubClassCode
static struct class_graph {
    const char *name;
    sipTypeDef **type;
    int yes, no;
} graph[] = {
    {sipName_QWinThumbnailToolBar, &sipType_QWinThumbnailToolBar, -1, 1},
    {sipName_QWinJumpList, &sipType_QWinJumpList, -1, 2},
    {sipName_QWinTaskbarButton, &sipType_QWinTaskbarButton, -1, 3},
    {sipName_QWinThumbnailToolButton, &sipType_QWinThumbnailToolButton, -1, 4},
    {sipName_QWinTaskbarProgress, &sipType_QWinTaskbarProgress, -1, -1},
};

int i = 0;

sipType = NULL;

do
{
    struct class_graph *cg = &graph[i];

    if (cg->name != NULL && sipCpp->inherits(cg->name))
    {
        sipType = *cg->type;
        i = cg->yes;
    }
    else
        i = cg->no;
}
while (i >= 0);
%End

public:
    explicit QWinJumpList(QObject *parent /TransferThis/ = 0);
    ~QWinJumpList();

    QString identifier() const;
    void setIdentifier(const QString &identifier);

    QWinJumpListCategory *recent() const;
    QWinJumpListCategory *frequent() const;
    QWinJumpListCategory *tasks() const;

    QList<QWinJumpListCategory *> categories() const;
    void addCategory(QWinJumpListCategory *category /Transfer/);
    QWinJumpListCategory *addCategory(const QString &title,
            const QList<QWinJumpListItem *> items /Transfer/ = QList<QWinJumpListItem *>());

public slots:
    void clear();

private:
    QWinJumpList(const QWinJumpList &);
};

%End
