""" various miscellaneous utilities """

def is_little_endian():
    """ am I little endian """
    import sys
    return sys.byteorder == 'little'

def exclusive(*args):
    count = sum([arg is not None for arg in args])
    return count == 1
