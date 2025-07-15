import sys

if sys.platform == "win32":
    from .tz.win import tzres as tzres, tzwin as tzwin, tzwinlocal as tzwinlocal
