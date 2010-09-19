def exclusive(*args):
    count = sum([arg is not None for arg in args])
    return count == 1

