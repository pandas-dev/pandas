// This is the SIP interface definition for QWinThumbnailToolBar.
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

class QWinThumbnailToolBar : public QObject
{
%TypeHeaderCode
#include <qwinthumbnailtoolbar.h>
%End

public:
    explicit QWinThumbnailToolBar(QObject *parent /TransferThis/ = 0);
    ~QWinThumbnailToolBar();

    void setWindow(QWindow *window);
    QWindow *window() const;

    void addButton(QWinThumbnailToolButton *button);
    void removeButton(QWinThumbnailToolButton *button);
    void setButtons(const QList<QWinThumbnailToolButton *> &buttons);
    QList<QWinThumbnailToolButton *> buttons() const;
    int count() const;

%If (Qt_5_4_0 -)
    bool iconicPixmapNotificationsEnabled() const;
    void setIconicPixmapNotificationsEnabled(bool enabled);

    QPixmap iconicThumbnailPixmap() const;
    QPixmap iconicLivePreviewPixmap() const;
%End

public slots:
    void clear();
%If (Qt_5_4_0 -)
    void setIconicThumbnailPixmap(const QPixmap &);
    void setIconicLivePreviewPixmap(const QPixmap &);
%End

signals:
%If (Qt_5_4_0 -)
    void iconicThumbnailPixmapRequested();
    void iconicLivePreviewPixmapRequested();
%End

private:
    QWinThumbnailToolBar(const QWinThumbnailToolBar &);
};

%End
