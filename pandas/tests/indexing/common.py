""" common utilities """


def _mklbl(prefix, n):
    return ["%s%s" % (prefix, i) for i in range(n)]
