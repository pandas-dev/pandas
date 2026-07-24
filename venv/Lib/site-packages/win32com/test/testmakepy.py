# Test makepy - try and run it over every OCX in the windows system directory.

import multiprocessing
import os
import shutil
import sys
import traceback

import pythoncom
import win32com.test.util
import winerror
from win32com.client import gencache, makepy, selecttlb


def TestBuildAll(verbose=1):
    num = 0
    tlbInfos = selecttlb.EnumTlbs()
    for info in tlbInfos:
        if verbose:
            print(f"{info.desc} ({info.dll})")
        try:
            makepy.GenerateFromTypeLibSpec(info)
            #          sys.stderr.write("Attr typeflags for coclass referenced object %s=%d (%d), typekind=%d\n" % (name, refAttr.wTypeFlags, refAttr.wTypeFlags & pythoncom.TYPEFLAG_FDUAL,refAttr.typekind))
            num += 1
        except pythoncom.com_error as details:
            # Ignore these 2 errors, as the are very common and can obscure
            # useful warnings.
            if details.hresult not in [
                winerror.TYPE_E_CANTLOADLIBRARY,
                winerror.TYPE_E_LIBNOTREGISTERED,
            ]:
                print("** COM error on", info.desc)
                print(details)
        except KeyboardInterrupt:
            print("Interrupted!")
            raise
        except:
            print("Failed:", info.desc)
            traceback.print_exc()
        if makepy.bForDemandDefault:
            # This only builds enums etc by default - build each
            # interface manually
            tinfo = (info.clsid, info.lcid, info.major, info.minor)
            mod = gencache.EnsureModule(info.clsid, info.lcid, info.major, info.minor)
            for name in mod.NamesToIIDMap:
                makepy.GenerateChildFromTypeLibSpec(name, tinfo)
    return num


def _TestEnsureModule(info):
    # This used to fail when called concurrently from multiple processes. See mhammond/pywin32#1923 .
    # The issue only happens when bForDemand is set as that creates a package instead of
    # a single module.
    tinfo = (info.clsid, info.lcid, int(info.major), int(info.minor))
    mod = gencache.EnsureModule(*tinfo, bForDemand=True)
    if makepy.bForDemandDefault:
        for name in mod.NamesToIIDMap:
            makepy.GenerateChildFromTypeLibSpec(name, tinfo)


def TestBuildConcurrent(verbose=1):
    # Pick any type library
    info = next(iter(selecttlb.EnumTlbs()))

    if verbose:
        print(f"{info.desc} ({info.dll})")

    # Call EnsureModule from multiple processes concurrently.
    nprocs = 16
    with multiprocessing.Pool(nprocs) as p:
        p.map(_TestEnsureModule, [info] * nprocs)

    return nprocs


def TestAll(verbose=0):
    gen_path = gencache.GetGeneratePath()
    if os.path.isdir(gen_path):
        shutil.rmtree(gen_path)

    nprocs = TestBuildConcurrent(verbose)
    print("Tested", nprocs, "concurrent processes")

    num = TestBuildAll(verbose)
    print("Generated and imported", num, "modules")

    win32com.test.util.CheckClean()


if __name__ == "__main__":
    TestAll("-q" not in sys.argv)
