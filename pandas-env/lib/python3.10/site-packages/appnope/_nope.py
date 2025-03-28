# -----------------------------------------------------------------------------
#  Copyright (C) 2013 Min RK
#
#  Distributed under the terms of the 2-clause BSD License.
# -----------------------------------------------------------------------------

from contextlib import contextmanager

import ctypes
import ctypes.util

objc = ctypes.cdll.LoadLibrary(ctypes.util.find_library("objc"))
_ = ctypes.cdll.LoadLibrary(ctypes.util.find_library("Foundation"))

void_p = ctypes.c_void_p
ull = ctypes.c_uint64

objc.objc_getClass.restype = void_p
objc.sel_registerName.restype = void_p
objc.objc_msgSend.restype = void_p
objc.objc_msgSend.argtypes = [void_p, void_p]

msg = objc.objc_msgSend


def _utf8(s):
    """ensure utf8 bytes"""
    if not isinstance(s, bytes):
        s = s.encode("utf8")
    return s


def n(name):
    """create a selector name (for methods)"""
    return objc.sel_registerName(_utf8(name))


def C(classname):
    """get an ObjC Class by name"""
    ret = objc.objc_getClass(_utf8(classname))
    assert ret is not None, "Couldn't find Class %s" % classname
    return ret


# constants from Foundation

NSActivityIdleDisplaySleepDisabled = 1 << 40
NSActivityIdleSystemSleepDisabled = 1 << 20
NSActivitySuddenTerminationDisabled = 1 << 14
NSActivityAutomaticTerminationDisabled = 1 << 15
NSActivityUserInitiated = 0x00FFFFFF | NSActivityIdleSystemSleepDisabled
NSActivityUserInitiatedAllowingIdleSystemSleep = (
    NSActivityUserInitiated & ~NSActivityIdleSystemSleepDisabled
)
NSActivityBackground = 0x000000FF
NSActivityLatencyCritical = 0xFF00000000


def beginActivityWithOptions(options, reason=""):
    """Wrapper for:

    [ [ NSProcessInfo processInfo]
        beginActivityWithOptions: (uint64)options
                          reason: (str)reason
    ]
    """
    NSProcessInfo = C("NSProcessInfo")
    NSString = C("NSString")

    objc.objc_msgSend.argtypes = [void_p, void_p, void_p]
    reason = msg(NSString, n("stringWithUTF8String:"), _utf8(reason))
    objc.objc_msgSend.argtypes = [void_p, void_p]
    info = msg(NSProcessInfo, n("processInfo"))
    objc.objc_msgSend.argtypes = [void_p, void_p, ull, void_p]
    activity = msg(
        info, n("beginActivityWithOptions:reason:"), ull(options), void_p(reason)
    )
    return activity


def endActivity(activity):
    """end a process activity assertion"""
    NSProcessInfo = C("NSProcessInfo")
    objc.objc_msgSend.argtypes = [void_p, void_p]
    info = msg(NSProcessInfo, n("processInfo"))
    objc.objc_msgSend.argtypes = [void_p, void_p, void_p]
    msg(info, n("endActivity:"), void_p(activity))


_theactivity = None


def nope():
    """disable App Nap by setting NSActivityUserInitiatedAllowingIdleSystemSleep"""
    global _theactivity
    _theactivity = beginActivityWithOptions(
        NSActivityUserInitiatedAllowingIdleSystemSleep, "Because Reasons"
    )


def nap():
    """end the caffeinated state started by `nope`"""
    global _theactivity
    if _theactivity is not None:
        endActivity(_theactivity)
        _theactivity = None


def napping_allowed():
    """is napping allowed?"""
    return _theactivity is None


@contextmanager
def nope_scope(
    options=NSActivityUserInitiatedAllowingIdleSystemSleep, reason="Because Reasons"
):
    """context manager for beginActivityWithOptions.

    Within this context, App Nap will be disabled.
    """
    activity = beginActivityWithOptions(options, reason)
    try:
        yield
    finally:
        endActivity(activity)


__all__ = [
    "NSActivityIdleDisplaySleepDisabled",
    "NSActivityIdleSystemSleepDisabled",
    "NSActivitySuddenTerminationDisabled",
    "NSActivityAutomaticTerminationDisabled",
    "NSActivityUserInitiated",
    "NSActivityUserInitiatedAllowingIdleSystemSleep",
    "NSActivityBackground",
    "NSActivityLatencyCritical",
    "beginActivityWithOptions",
    "endActivity",
    "nope",
    "nap",
    "napping_allowed",
    "nope_scope",
]
