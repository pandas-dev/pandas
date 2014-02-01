#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
import re
import os
import time
import json
import base64

from iron_cache import *

key='KEY'+os.environ.get('JOB_NAME','')+"."

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
