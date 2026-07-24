// This is the SIP interface definition for QtWin.
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

typedef struct HBITMAP__ *HBITMAP;
typedef struct HDC__ *HDC;
typedef struct HICON__ *HICON;
typedef long HRESULT;
typedef struct HRGN__ *HRGN;


namespace QtWin
{
%TypeHeaderCode
#include <qwinfunctions.h>
%End

    enum HBitmapFormat
    {
        HBitmapNoAlpha,
        HBitmapPremultipliedAlpha,
        HBitmapAlpha
    };

    enum WindowFlip3DPolicy
    {
        FlipDefault,
        FlipExcludeBelow,
        FlipExcludeAbove
    };

    HBITMAP createMask(const QBitmap &bitmap);
    HBITMAP toHBITMAP(const QPixmap &p,
            QtWin::HBitmapFormat format = QtWin::HBitmapNoAlpha);
    QPixmap fromHBITMAP(HBITMAP bitmap,
            QtWin::HBitmapFormat format = QtWin::HBitmapNoAlpha);
    HICON toHICON(const QPixmap &p);
    QImage imageFromHBITMAP(HDC hdc, HBITMAP bitmap, int width, int height);
    QPixmap fromHICON(HICON icon);
    HRGN toHRGN(const QRegion &region);
    QRegion fromHRGN(HRGN hrgn);

    QString stringFromHresult(HRESULT hresult);
    QString errorStringFromHresult(HRESULT hresult);

    QColor colorizationColor(bool *opaqueBlend = 0);
    QColor realColorizationColor();

    void setWindowExcludedFromPeek(QWindow *window, bool exclude);
    bool isWindowExcludedFromPeek(QWindow *window);
    void setWindowDisallowPeek(QWindow *window, bool disallow);
    bool isWindowPeekDisallowed(QWindow *window);
    void setWindowFlip3DPolicy(QWindow *window,
            QtWin::WindowFlip3DPolicy policy);
    QtWin::WindowFlip3DPolicy windowFlip3DPolicy(QWindow *);

    void extendFrameIntoClientArea(QWindow *window, int left, int top,
            int right, int bottom);
    void extendFrameIntoClientArea(QWindow *window, const QMargins &margins);
    void resetExtendedFrame(QWindow *window);

    void enableBlurBehindWindow(QWindow *window, const QRegion &region);
    void enableBlurBehindWindow(QWindow *window);
    void disableBlurBehindWindow(QWindow *window);

    bool isCompositionEnabled();
    void setCompositionEnabled(bool enabled);
    bool isCompositionOpaque();

    void setCurrentProcessExplicitAppUserModelID(const QString &id);

    void markFullscreenWindow(QWindow *, bool fullscreen = true);

    void taskbarActivateTab(QWindow *);
    void taskbarActivateTabAlt(QWindow *);
    void taskbarAddTab(QWindow *);
    void taskbarDeleteTab(QWindow *);

    void setWindowExcludedFromPeek(QWidget *window, bool exclude);
    bool isWindowExcludedFromPeek(QWidget *window);
    void setWindowDisallowPeek(QWidget *window, bool disallow);
    bool isWindowPeekDisallowed(QWidget *window);
    void setWindowFlip3DPolicy(QWidget *window,
            QtWin::WindowFlip3DPolicy policy);
    QtWin::WindowFlip3DPolicy windowFlip3DPolicy(QWidget *window);

    void extendFrameIntoClientArea(QWidget *window, const QMargins &margins);
    void extendFrameIntoClientArea(QWidget *window, int left, int top,
            int right, int bottom);
    void resetExtendedFrame(QWidget *window);

    void enableBlurBehindWindow(QWidget *window, const QRegion &region);
    void enableBlurBehindWindow(QWidget *window);
    void disableBlurBehindWindow(QWidget *window);

    void markFullscreenWindow(QWidget *window, bool fullscreen = true);

    void taskbarActivateTab(QWidget *window);
    void taskbarActivateTabAlt(QWidget *window);
    void taskbarAddTab(QWidget *window);
    void taskbarDeleteTab(QWidget *window);
};

%End
