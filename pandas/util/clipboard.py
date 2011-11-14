"""
Taken from the IPython project http://ipython.org

Used under the terms of the BSD license
"""

import subprocess
import sys

def clipboard_get():
    """ Get text from the clipboard.
    """
    if sys.platform == 'win32':
        try:
            return win32_clipboard_get()
        except Exception:
            pass
    elif sys.platform == 'darwin':
        try:
            return osx_clipboard_get()
        except Exception:
            pass
    return tkinter_clipboard_get()

def win32_clipboard_get():
    """ Get the current clipboard's text on Windows.

    Requires Mark Hammond's pywin32 extensions.
    """
    try:
        import win32clipboard
    except ImportError:
        message = ("Getting text from the clipboard requires the pywin32 "
            "extensions: http://sourceforge.net/projects/pywin32/")
        raise Exception(message)
    win32clipboard.OpenClipboard()
    text = win32clipboard.GetClipboardData(win32clipboard.CF_TEXT)
    # FIXME: convert \r\n to \n?
    win32clipboard.CloseClipboard()
    return text

def osx_clipboard_get():
    """ Get the clipboard's text on OS X.
    """
    p = subprocess.Popen(['pbpaste', '-Prefer', 'ascii'],
        stdout=subprocess.PIPE)
    text, stderr = p.communicate()
    # Text comes in with old Mac \r line endings. Change them to \n.
    text = text.replace('\r', '\n')
    return text

def tkinter_clipboard_get():
    """ Get the clipboard's text using Tkinter.

    This is the default on systems that are not Windows or OS X. It may
    interfere with other UI toolkits and should be replaced with an
    implementation that uses that toolkit.
    """
    try:
        import Tkinter
    except ImportError:
        message = ("Getting text from the clipboard on this platform "
            "requires Tkinter.")
        raise Exception(message)
    root = Tkinter.Tk()
    root.withdraw()
    text = root.clipboard_get()
    root.destroy()
    return text
