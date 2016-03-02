# Pyperclip v1.5.15
# A cross-platform clipboard module for Python.
# By Al Sweigart al@inventwithpython.com

# Usage:
#   import pyperclip
#   pyperclip.copy('The text to be copied to the clipboard.')
#   spam = pyperclip.paste()

# On Windows, no additional modules are needed.
# On Mac, this module makes use of the pbcopy and pbpaste commands, which
# should come with the os.
# On Linux, this module makes use of the xclip or xsel commands, which should
# come with the os. Otherwise run "sudo apt-get install xclip" or
# "sudo apt-get install xsel"
# Otherwise on Linux, you will need the gtk or PyQt4 modules installed.
# The gtk module is not available for Python 3, and this module does not work
# with PyGObject yet.


# Copyright (c) 2015, Albert Sweigart
# All rights reserved.
#
# BSD-style license:
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
#     * Redistributions of source code must retain the above copyright
#       notice, this list of conditions and the following disclaimer.
#     * Redistributions in binary form must reproduce the above copyright
#       notice, this list of conditions and the following disclaimer in the
#       documentation and/or other materials provided with the distribution.
#     * Neither the name of the pyperclip nor the
#       names of its contributors may be used to endorse or promote products
#       derived from this software without specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY Albert Sweigart "AS IS" AND ANY
# EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
# WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
# DISCLAIMED. IN NO EVENT SHALL Albert Sweigart BE LIABLE FOR ANY
# DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
# (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
# LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
# ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
# (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
# SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

# flake8: noqa

import platform
import os
from subprocess import call, Popen, PIPE

PY2 = '2' == platform.python_version_tuple()[0]
text_type = unicode if PY2 else str


class NoClipboardProgramError(OSError):
    pass


def _pasteWindows():
    CF_UNICODETEXT = 13
    d = ctypes.windll
    d.user32.OpenClipboard(0)
    handle = d.user32.GetClipboardData(CF_UNICODETEXT)
    data = ctypes.c_wchar_p(handle).value
    d.user32.CloseClipboard()
    return data


def _copyWindows(text):
    GMEM_DDESHARE = 0x2000
    CF_UNICODETEXT = 13
    d = ctypes.windll  # cdll expects 4 more bytes in user32.OpenClipboard(0)
    if not isinstance(text, text_type):
        text = text.decode('mbcs')

    d.user32.OpenClipboard(0)

    d.user32.EmptyClipboard()
    hCd = d.kernel32.GlobalAlloc(GMEM_DDESHARE,
                                 len(text.encode('utf-16-le')) + 2)
    pchData = d.kernel32.GlobalLock(hCd)
    ctypes.cdll.msvcrt.wcscpy(ctypes.c_wchar_p(pchData), text)
    d.kernel32.GlobalUnlock(hCd)
    d.user32.SetClipboardData(CF_UNICODETEXT, hCd)
    d.user32.CloseClipboard()


def _pasteCygwin():
    CF_UNICODETEXT = 13
    d = ctypes.cdll
    d.user32.OpenClipboard(0)
    handle = d.user32.GetClipboardData(CF_UNICODETEXT)
    data = ctypes.c_wchar_p(handle).value
    d.user32.CloseClipboard()
    return data


def _copyCygwin(text):
    GMEM_DDESHARE = 0x2000
    CF_UNICODETEXT = 13
    d = ctypes.cdll
    if not isinstance(text, text_type):
        text = text.decode('mbcs')
    d.user32.OpenClipboard(0)
    d.user32.EmptyClipboard()
    hCd = d.kernel32.GlobalAlloc(GMEM_DDESHARE,
                                 len(text.encode('utf-16-le')) + 2)
    pchData = d.kernel32.GlobalLock(hCd)
    ctypes.cdll.msvcrt.wcscpy(ctypes.c_wchar_p(pchData), text)
    d.kernel32.GlobalUnlock(hCd)
    d.user32.SetClipboardData(CF_UNICODETEXT, hCd)
    d.user32.CloseClipboard()


def _copyOSX(text):
    p = Popen(['pbcopy', 'w'], stdin=PIPE, close_fds=True)
    p.communicate(input=text.encode('utf-8'))


def _pasteOSX():
    p = Popen(['pbpaste', 'r'], stdout=PIPE, close_fds=True)
    stdout, stderr = p.communicate()
    return stdout.decode('utf-8')


def _pasteGtk():
    return gtk.Clipboard().wait_for_text()


def _copyGtk(text):
    global cb
    cb = gtk.Clipboard()
    cb.set_text(text)
    cb.store()


def _pasteQt():
    return str(cb.text())


def _copyQt(text):
    cb.setText(text)


def _copyXclip(text):
    p = Popen(['xclip', '-selection', 'c'], stdin=PIPE, close_fds=True)
    p.communicate(input=text.encode('utf-8'))


def _pasteXclip():
    p = Popen(['xclip', '-selection', 'c', '-o'], stdout=PIPE, close_fds=True)
    stdout, stderr = p.communicate()
    return stdout.decode('utf-8')


def _copyXsel(text):
    p = Popen(['xsel', '-b', '-i'], stdin=PIPE, close_fds=True)
    p.communicate(input=text.encode('utf-8'))


def _pasteXsel():
    p = Popen(['xsel', '-b', '-o'], stdout=PIPE, close_fds=True)
    stdout, stderr = p.communicate()
    return stdout.decode('utf-8')


def _copyKlipper(text):
    p = Popen(['qdbus', 'org.kde.klipper', '/klipper',
               'setClipboardContents', text.encode('utf-8')],
              stdin=PIPE, close_fds=True)
    p.communicate(input=None)


def _pasteKlipper():
    p = Popen(['qdbus', 'org.kde.klipper', '/klipper',
               'getClipboardContents'], stdout=PIPE, close_fds=True)
    stdout, stderr = p.communicate()
    return stdout.decode('utf-8')


# Determine the OS/platform and set the copy() and paste() functions
# accordingly.
if 'cygwin' in platform.system().lower():
    _functions = 'Cygwin'  # for debugging
    import ctypes
    paste = _pasteCygwin
    copy = _copyCygwin
elif os.name == 'nt' or platform.system() == 'Windows':
    _functions = 'Windows'  # for debugging
    import ctypes
    paste = _pasteWindows
    copy = _copyWindows
elif os.name == 'mac' or platform.system() == 'Darwin':
    _functions = 'OS X pbcopy/pbpaste'  # for debugging
    paste = _pasteOSX
    copy = _copyOSX
elif os.name == 'posix' or platform.system() == 'Linux':
    # Determine which command/module is installed, if any.
    xclipExists = call(['which', 'xclip'],
                       stdout=PIPE, stderr=PIPE) == 0

    xselExists = call(['which', 'xsel'],
                      stdout=PIPE, stderr=PIPE) == 0

    xklipperExists = (
        call(['which', 'klipper'], stdout=PIPE, stderr=PIPE) == 0 and
        call(['which', 'qdbus'], stdout=PIPE, stderr=PIPE) == 0
    )

    gtkInstalled = False
    try:
        # Check it gtk is installed.
        import gtk
        gtkInstalled = True
    except ImportError:
        pass

    if not gtkInstalled:
        # Check for either PyQt4 or PySide
        qtBindingInstalled = True
        try:
            from PyQt4 import QtGui
        except ImportError:
            try:
                from PySide import QtGui
            except ImportError:
                qtBindingInstalled = False

    # Set one of the copy & paste functions.
    if xclipExists:
        _functions = 'xclip command'  # for debugging
        paste = _pasteXclip
        copy = _copyXclip
    elif xklipperExists:
        _functions = '(KDE Klipper) - qdbus (external)'  # for debugging
        paste = _pasteKlipper
        copy = _copyKlipper
    elif gtkInstalled:
        _functions = 'gtk module'  # for debugging
        paste = _pasteGtk
        copy = _copyGtk
    elif qtBindingInstalled:
        _functions = 'PyQt4 module'  # for debugging
        app = QtGui.QApplication([])
        cb = QtGui.QApplication.clipboard()
        paste = _pasteQt
        copy = _copyQt
    elif xselExists:
        # TODO: xsel doesn't seem to work on Raspberry Pi (my test Linux
        # environment). Putting this as the last method tried.
        _functions = 'xsel command'  # for debugging
        paste = _pasteXsel
        copy = _copyXsel
    else:
        raise NoClipboardProgramError('Pyperclip requires the gtk, PyQt4, or '
                                      'PySide module installed, or either the '
                                      'xclip or xsel command.')
else:
    raise RuntimeError('pyperclip does not support your system.')

# pandas aliases
clipboard_get = paste
clipboard_set = copy
