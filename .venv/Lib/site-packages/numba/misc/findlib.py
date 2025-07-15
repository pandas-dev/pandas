import sys
import os
import re


def get_lib_dirs():
    """
    Anaconda specific
    """
    if sys.platform == 'win32':
        # on windows, historically `DLLs` has been used for CUDA libraries,
        # since approximately CUDA 9.2, `Library\bin` has been used.
        dirnames = ['DLLs', os.path.join('Library', 'bin')]
    else:
        dirnames = ['lib', ]
    libdirs = [os.path.join(sys.prefix, x) for x in dirnames]
    return libdirs


DLLNAMEMAP = {
    'linux': r'lib%(name)s\.so\.%(ver)s$',
    'linux2': r'lib%(name)s\.so\.%(ver)s$',
    'linux-static': r'lib%(name)s\.a$',
    'darwin': r'lib%(name)s\.%(ver)s\.dylib$',
    'win32': r'%(name)s%(ver)s\.dll$',
    'win32-static': r'%(name)s\.lib$',
    'bsd': r'lib%(name)s\.so\.%(ver)s$',
}

RE_VER = r'[0-9]*([_\.][0-9]+)*'


def find_lib(libname, libdir=None, platform=None, static=False):
    platform = platform or sys.platform
    platform = 'bsd' if 'bsd' in platform else platform
    if static:
        platform = f"{platform}-static"
    if platform not in DLLNAMEMAP:
        # Return empty list if platform name is undefined.
        # Not all platforms define their static library paths.
        return []
    pat = DLLNAMEMAP[platform] % {"name": libname, "ver": RE_VER}
    regex = re.compile(pat)
    return find_file(regex, libdir)


def find_file(pat, libdir=None):
    if libdir is None:
        libdirs = get_lib_dirs()
    elif isinstance(libdir, str):
        libdirs = [libdir,]
    else:
        libdirs = list(libdir)
    files = []
    for ldir in libdirs:
        try:
            entries = os.listdir(ldir)
        except FileNotFoundError:
            continue
        candidates = [os.path.join(ldir, ent)
                      for ent in entries if pat.match(ent)]
        files.extend([c for c in candidates if os.path.isfile(c)])
    return files
