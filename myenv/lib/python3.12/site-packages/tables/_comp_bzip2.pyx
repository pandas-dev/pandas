import sys
from libc.stdlib cimport free


cdef extern from "H5Zbzip2.h":
    int register_bzip2(char **, char **)


def register_():
    cdef char *version
    cdef char *date

    if not register_bzip2(&version, &date):
        return None

    compinfo = (version, date)
    free(version)
    free(date)
    return compinfo[0].decode('ascii'), compinfo[1].decode('ascii')
