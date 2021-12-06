# OrderedDict is kind of slow, so for most of our uses in Python 3.6
# and later we'd rather just use dict

import sys

if sys.version_info < (3, 6):
    from collections import OrderedDict as OrderedDict
else:
    OrderedDict = dict
