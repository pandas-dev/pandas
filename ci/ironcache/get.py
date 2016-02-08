#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
import re
import os
import time
import json
import base64
from hashlib import sha1
from iron_cache import *
import traceback as tb

key='KEY.%s.%s' %(os.environ.get('TRAVIS_REPO_SLUG','unk'),
                   os.environ.get('JOB_NAME','unk'))
print(key)

if sys.version_info[0] > 2:
    key = bytes(key,encoding='utf8')

key = sha1(key).hexdigest()[:8]+'.'

b = b''
cache = IronCache()
for i in range(20):
    print("getting %s" % key+str(i))
    try:
        item = cache.get(cache="travis", key=key+str(i))
        v = item.value
        if sys.version_info[0] > 2:
            v = bytes(v,encoding='utf8')
        b += bytes(base64.b64decode(v))
    except Exception as e:
        try:
            print(tb.format_exc(e))
        except:
            print("exception during exception, oh my")
        break

with open(os.path.join(os.environ.get('HOME',''),"ccache.7z"),'wb') as f:
    f.write(b)
