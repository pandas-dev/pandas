import numpy as np


def _ensure_decoded(s):
    """ if we have bytes, decode them to unicode """
    if isinstance(s, (np.bytes_, bytes)):
        s = s.decode('UTF-8')
    return s


class NameResolutionError(NameError):
    pass
