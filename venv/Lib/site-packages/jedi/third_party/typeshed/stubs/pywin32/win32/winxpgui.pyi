# The `winxpgui` module is obsolete and has been completely replaced
# by `win32gui` and `win32console.GetConsoleWindow`. Use those instead.

from win32console import GetConsoleWindow as GetConsoleWindow
from win32gui import *
