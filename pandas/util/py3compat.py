import sys

PY3 = (sys.version_info[0] >= 3)

if PY3:
    def isidentifier(s):
        return s.isidentifier()

    def str_to_bytes(s, encoding='ascii'):
        return s.encode(encoding)

    def bytes_to_str(b, encoding='utf-8'):
        return b.decode(encoding)

    lzip = lambda *args: list(zip(*args))
else:
    # Python 2
    import re
    _name_re = re.compile(r"[a-zA-Z_][a-zA-Z0-9_]*$")

    def isidentifier(s, dotted=False):
        return bool(_name_re.match(s))

    def str_to_bytes(s, encoding='ascii'):
        return s

    def bytes_to_str(b, encoding='ascii'):
        return b

    lzip = zip

try:
    from cStringIO import StringIO
except:
    from io import StringIO

try:
    from io import BytesIO
except:
    from cStringIO import StringIO as BytesIO
