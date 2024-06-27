// This is the SIP specification of the QSysInfo class.
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


class QSysInfo
{
%TypeHeaderCode
#include <qsysinfo.h>
%End

public:
    enum Sizes
    {
        WordSize,
    };

    enum Endian
    {
        BigEndian,
        LittleEndian,
        ByteOrder,
    };

%If (WS_MACX)
    enum MacVersion
    {
        MV_Unknown,
        MV_9,
        MV_10_0,
        MV_10_1,
        MV_10_2,
        MV_10_3,
        MV_10_4,
        MV_10_5,
        MV_10_6,
        MV_10_7,
        MV_10_8,
%If (Qt_5_1_0 -)
        MV_10_9,
%End
%If (Qt_5_4_0 -)
        MV_10_10,
%End
%If (Qt_5_5_0 -)
        MV_10_11,
%End
%If (Qt_5_8_0 -)
        MV_10_12,
%End
        MV_CHEETAH,
        MV_PUMA,
        MV_JAGUAR,
        MV_PANTHER,
        MV_TIGER,
        MV_LEOPARD,
        MV_SNOWLEOPARD,
        MV_LION,
        MV_MOUNTAINLION,
%If (Qt_5_1_0 -)
        MV_MAVERICKS,
%End
%If (Qt_5_4_0 -)
        MV_YOSEMITE,
%End
%If (Qt_5_5_0 -)
        MV_ELCAPITAN,
%End
%If (Qt_5_8_0 -)
        MV_SIERRA,
%End

%If (Qt_5_2_0 -)
        MV_IOS,
        MV_IOS_4_3,
        MV_IOS_5_0,
        MV_IOS_5_1,
        MV_IOS_6_0,
        MV_IOS_6_1,
        MV_IOS_7_0,
        MV_IOS_7_1,
%End
%If (Qt_5_4_0 -)
        MV_IOS_8_0,
%End
%If (Qt_5_5_0 -)
        MV_IOS_8_1,
        MV_IOS_8_2,
        MV_IOS_8_3,
        MV_IOS_8_4,
        MV_IOS_9_0,
%End
%If (Qt_5_8_0 -)
        MV_IOS_9_1,
        MV_IOS_9_2,
        MV_IOS_9_3,
        MV_IOS_10_0,
%End

%If (Qt_5_8_0 -)
        MV_TVOS,
        MV_TVOS_9_0,
        MV_TVOS_9_1,
        MV_TVOS_9_2,
        MV_TVOS_10_0,
%End

%If (Qt_5_8_0 -)
        MV_WATCHOS,
        MV_WATCHOS_2_0,
        MV_WATCHOS_2_1,
        MV_WATCHOS_2_2,
        MV_WATCHOS_3_0,
%End
    };

    static const QSysInfo::MacVersion MacintoshVersion;
    static QSysInfo::MacVersion macVersion();
%End

%If (WS_WIN)
    enum WinVersion {
        WV_32s,
        WV_95,
        WV_98,
        WV_Me,
        WV_DOS_based,

        WV_NT,
        WV_2000,
        WV_XP,
        WV_2003,
        WV_VISTA,
        WV_WINDOWS7,
        WV_WINDOWS8,
%If (Qt_5_2_0 -)
        WV_WINDOWS8_1,
%End
%If (Qt_5_5_0 -)
        WV_WINDOWS10,
%End
        WV_NT_based,

        WV_4_0,
        WV_5_0,
        WV_5_1,
        WV_5_2,
        WV_6_0,
        WV_6_1,
        WV_6_2,
%If (Qt_5_2_0 -)
        WV_6_3,
%End
%If (Qt_5_5_0 -)
        WV_10_0,
%End

        WV_CE,
        WV_CENET,
        WV_CE_5,
        WV_CE_6,
        WV_CE_based
    };

    static const QSysInfo::WinVersion WindowsVersion;
    static QSysInfo::WinVersion windowsVersion();
%End

%If (Qt_5_4_0 -)
    static QString buildAbi();
    static QString buildCpuArchitecture();
    static QString currentCpuArchitecture();
    static QString kernelType();
    static QString kernelVersion();
    static QString prettyProductName();
    static QString productType();
    static QString productVersion();
%End

%If (Qt_5_6_0 -)
    static QString machineHostName();
%End
};
