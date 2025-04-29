"""
Replacement module for what html5lib uses six for.
"""

import http.client
import operator
import urllib


PY3 = True
binary_type = bytes
string_types = (str,)
text_type = str
unichr = chr
viewkeys = operator.methodcaller("keys")

http_client = http.client
urllib = urllib
urllib_parse = urllib.parse
