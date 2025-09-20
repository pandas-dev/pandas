// This is the SIP interface definition for QX11Info.
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


%If (Qt_5_1_0 -)

class Display;
class xcb_connection_t;


class QX11Info
{
%TypeHeaderCode
#include <qx11info_x11.h>
%End

public:
%If (Qt_5_2_0 -)
    static bool isPlatformX11();
%End

    static int appDpiX(int screen = -1);
    static int appDpiY(int screen = -1);

    static unsigned long appRootWindow(int screen = -1);
    static int appScreen();

    static unsigned long appTime();
    static unsigned long appUserTime();

    static void setAppTime(unsigned long time);
    static void setAppUserTime(unsigned long time);

%If (Qt_5_2_0 -)
    static unsigned long getTimestamp();
%End

%If (Qt_5_4_0 -)
    static QByteArray nextStartupId();
    static void setNextStartupId(const QByteArray &id);
%End

    static Display *display();
    static xcb_connection_t *connection();

private:
    QX11Info();
};

%End
