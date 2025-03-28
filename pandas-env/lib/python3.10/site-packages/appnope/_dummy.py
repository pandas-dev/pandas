# -----------------------------------------------------------------------------
#  Copyright (C) 2013 Min RK
#
#  Distributed under the terms of the 2-clause BSD License.
# -----------------------------------------------------------------------------

from contextlib import contextmanager


def beginActivityWithOptions(options, reason=""):
    return


def endActivity(activity):
    return


def nope():
    return


def nap():
    return


@contextmanager
def nope_scope(options=0, reason="Because Reasons"):
    yield


def napping_allowed():
    return True
