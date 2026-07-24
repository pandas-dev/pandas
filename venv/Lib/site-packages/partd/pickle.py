"""
get/put functions that consume/produce Python lists using Pickle to serialize
"""
import pickle

from .encode import Encode
from functools import partial

def concat(lists):
    return sum(lists, [])

Pickle = partial(Encode,
                 partial(pickle.dumps, protocol=pickle.HIGHEST_PROTOCOL),
                 pickle.loads,
                 concat)
