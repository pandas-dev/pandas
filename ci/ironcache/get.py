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

key='KEY.%s.%s' %(os.environ.get('TRAVIS_REPO_SLUG','unk'),
                   os.environ.get('JOB_NAME','unk'))
if sys.version_info[0] > 2:
    key = sha1(bytes(key,encoding='utf8')).hexdigest()[:8]+'.'
else:
    key = sha1(key).hexdigest()[:8]+'.'

b = b''
cache = IronCache()
for i in range(20):
    print("getting %s" % key+str(i))
    try:
        item = cache.get(cache="travis", key=key+str(i))
        b += bytes(base64.b64decode(item.value))
    except:
        break

with open(os.path.join(os.environ.get('HOME',''),"ccache.7z"),'wb') as f:
    f.write(b)
