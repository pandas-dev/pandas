import sys
if sys.platform == 'win32':
    from pandas._libs.src.tzlocal.win32 import get_localzone, reload_localzone
else:
    from pandas._libs.src.tzlocal.unix import get_localzone, reload_localzone
