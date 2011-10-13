import sys

PY3 = (sys.version_info[0] >= 3)

if PY3:
    def isidentifier(s):
        return s.isidentifier()

else:
    # Python 2
    import re
    _name_re = re.compile(r"[a-zA-Z_][a-zA-Z0-9_]*$")
    def isidentifier(s, dotted=False):
        return bool(_name_re.match(s))
