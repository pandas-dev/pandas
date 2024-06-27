#############################################################################
##
## Copyright (c) 2023 Riverbank Computing Limited <info@riverbankcomputing.com>
## 
## This file is part of PyQt5.
## 
## This file may be used under the terms of the GNU General Public License
## version 3.0 as published by the Free Software Foundation and appearing in
## the file LICENSE included in the packaging of this file.  Please review the
## following information to ensure the GNU General Public License version 3.0
## requirements will be met: http://www.gnu.org/copyleft/gpl.html.
## 
## If you do not wish to use this file under the terms of the GPL version 3.0
## then you may purchase a commercial license.  For more information contact
## info@riverbankcomputing.com.
## 
## This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
## WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
##
#############################################################################


import re


def as_string(obj):
    if isinstance(obj, basestring):
        return '"' + _escape(obj.encode('UTF-8')) + '"'

    return str(obj)


_esc_regex = re.compile(r"(\"|\'|\\)")

def _escape(text):
    # This escapes any escaped single or double quote or backslash.
    x = _esc_regex.sub(r"\\\1", text)

    # This replaces any '\n' with an escaped version and a real line break.
    return re.sub(r'\n', r'\\n"\n"', x)
