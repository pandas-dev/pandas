// This is the SIP interface definition for QWinTaskbarProgress.
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

class QWinTaskbarProgress : public QObject
{
%TypeHeaderCode
#include <qwintaskbarprogress.h>
%End

public:
    explicit QWinTaskbarProgress(QObject *parent /TransferThis/= 0);
    ~QWinTaskbarProgress();

    int value() const;
    int minimum() const;
    int maximum() const;
    bool isVisible() const;
    bool isPaused() const;
    bool isStopped() const;

public slots:
    void setValue(int value);
    void setMinimum(int minimum);
    void setMaximum(int maximum);
    void setRange(int minimum, int maximum);
    void reset();
    void show();
    void hide();
    void setVisible(bool visible);
    void pause();
    void resume();
    void setPaused(bool paused);
    void stop();

signals:
    void valueChanged(int value);
    void minimumChanged(int minimum);
    void maximumChanged(int maximum);
    void visibilityChanged(bool visible);

private:
    QWinTaskbarProgress(const QWinTaskbarProgress &);
};

%End
