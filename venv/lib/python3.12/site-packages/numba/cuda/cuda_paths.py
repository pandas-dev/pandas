import sys
import re
import os
from collections import namedtuple

from numba.core.config import IS_WIN32
from numba.misc.findlib import find_lib, find_file


_env_path_tuple = namedtuple('_env_path_tuple', ['by', 'info'])


def _find_valid_path(options):
    """Find valid path from *options*, which is a list of 2-tuple of
    (name, path).  Return first pair where *path* is not None.
    If no valid path is found, return ('<unknown>', None)
    """
    for by, data in options:
        if data is not None:
            return by, data
    else:
        return '<unknown>', None


def _get_libdevice_path_decision():
    options = [
        ('Conda environment', get_conda_ctk()),
        ('Conda environment (NVIDIA package)', get_nvidia_libdevice_ctk()),
        ('CUDA_HOME', get_cuda_home('nvvm', 'libdevice')),
        ('System', get_system_ctk('nvvm', 'libdevice')),
        ('Debian package', get_debian_pkg_libdevice()),
    ]
    by, libdir = _find_valid_path(options)
    return by, libdir


def _nvvm_lib_dir():
    if IS_WIN32:
        return 'nvvm', 'bin'
    else:
        return 'nvvm', 'lib64'


def _get_nvvm_path_decision():
    options = [
        ('Conda environment', get_conda_ctk()),
        ('Conda environment (NVIDIA package)', get_nvidia_nvvm_ctk()),
        ('CUDA_HOME', get_cuda_home(*_nvvm_lib_dir())),
        ('System', get_system_ctk(*_nvvm_lib_dir())),
    ]
    by, path = _find_valid_path(options)
    return by, path


def _get_libdevice_paths():
    by, libdir = _get_libdevice_path_decision()
    # Search for pattern
    pat = r'libdevice(\.\d+)*\.bc$'
    candidates = find_file(re.compile(pat), libdir)
    # Keep only the max (most recent version) of the bitcode files.
    out = max(candidates, default=None)
    return _env_path_tuple(by, out)


def _cudalib_path():
    if IS_WIN32:
        return 'bin'
    else:
        return 'lib64'


def _cuda_home_static_cudalib_path():
    if IS_WIN32:
        return ('lib', 'x64')
    else:
        return ('lib64',)


def _get_cudalib_dir_path_decision():
    options = [
        ('Conda environment', get_conda_ctk()),
        ('Conda environment (NVIDIA package)', get_nvidia_cudalib_ctk()),
        ('CUDA_HOME', get_cuda_home(_cudalib_path())),
        ('System', get_system_ctk(_cudalib_path())),
    ]
    by, libdir = _find_valid_path(options)
    return by, libdir


def _get_static_cudalib_dir_path_decision():
    options = [
        ('Conda environment', get_conda_ctk()),
        ('Conda environment (NVIDIA package)', get_nvidia_static_cudalib_ctk()),
        ('CUDA_HOME', get_cuda_home(*_cuda_home_static_cudalib_path())),
        ('System', get_system_ctk(_cudalib_path())),
    ]
    by, libdir = _find_valid_path(options)
    return by, libdir


def _get_cudalib_dir():
    by, libdir = _get_cudalib_dir_path_decision()
    return _env_path_tuple(by, libdir)


def _get_static_cudalib_dir():
    by, libdir = _get_static_cudalib_dir_path_decision()
    return _env_path_tuple(by, libdir)


def get_system_ctk(*subdirs):
    """Return path to system-wide cudatoolkit; or, None if it doesn't exist.
    """
    # Linux?
    if sys.platform.startswith('linux'):
        # Is cuda alias to /usr/local/cuda?
        # We are intentionally not getting versioned cuda installation.
        base = '/usr/local/cuda'
        if os.path.exists(base):
            return os.path.join(base, *subdirs)


def get_conda_ctk():
    """Return path to directory containing the shared libraries of cudatoolkit.
    """
    is_conda_env = os.path.exists(os.path.join(sys.prefix, 'conda-meta'))
    if not is_conda_env:
        return
    # Assume the existence of NVVM to imply cudatoolkit installed
    paths = find_lib('nvvm')
    if not paths:
        return
    # Use the directory name of the max path
    return os.path.dirname(max(paths))


def get_nvidia_nvvm_ctk():
    """Return path to directory containing the NVVM shared library.
    """
    is_conda_env = os.path.exists(os.path.join(sys.prefix, 'conda-meta'))
    if not is_conda_env:
        return

    # Assume the existence of NVVM in the conda env implies that a CUDA toolkit
    # conda package is installed.

    # First, try the location used on Linux and the Windows 11.x packages
    libdir = os.path.join(sys.prefix, 'nvvm', _cudalib_path())
    if not os.path.exists(libdir) or not os.path.isdir(libdir):
        # If that fails, try the location used for Windows 12.x packages
        libdir = os.path.join(sys.prefix, 'Library', 'nvvm', _cudalib_path())
        if not os.path.exists(libdir) or not os.path.isdir(libdir):
            # If that doesn't exist either, assume we don't have the NVIDIA
            # conda package
            return

    paths = find_lib('nvvm', libdir=libdir)
    if not paths:
        return
    # Use the directory name of the max path
    return os.path.dirname(max(paths))


def get_nvidia_libdevice_ctk():
    """Return path to directory containing the libdevice library.
    """
    nvvm_ctk = get_nvidia_nvvm_ctk()
    if not nvvm_ctk:
        return
    nvvm_dir = os.path.dirname(nvvm_ctk)
    return os.path.join(nvvm_dir, 'libdevice')


def get_nvidia_cudalib_ctk():
    """Return path to directory containing the shared libraries of cudatoolkit.
    """
    nvvm_ctk = get_nvidia_nvvm_ctk()
    if not nvvm_ctk:
        return
    env_dir = os.path.dirname(os.path.dirname(nvvm_ctk))
    subdir = 'bin' if IS_WIN32 else 'lib'
    return os.path.join(env_dir, subdir)


def get_nvidia_static_cudalib_ctk():
    """Return path to directory containing the static libraries of cudatoolkit.
    """
    nvvm_ctk = get_nvidia_nvvm_ctk()
    if not nvvm_ctk:
        return

    if IS_WIN32 and ("Library" not in nvvm_ctk):
        # Location specific to CUDA 11.x packages on Windows
        dirs = ('Lib', 'x64')
    else:
        # Linux, or Windows with CUDA 12.x packages
        dirs = ('lib',)

    env_dir = os.path.dirname(os.path.dirname(nvvm_ctk))
    return os.path.join(env_dir, *dirs)


def get_cuda_home(*subdirs):
    """Get paths of CUDA_HOME.
    If *subdirs* are the subdirectory name to be appended in the resulting
    path.
    """
    cuda_home = os.environ.get('CUDA_HOME')
    if cuda_home is None:
        # Try Windows CUDA installation without Anaconda
        cuda_home = os.environ.get('CUDA_PATH')
    if cuda_home is not None:
        return os.path.join(cuda_home, *subdirs)


def _get_nvvm_path():
    by, path = _get_nvvm_path_decision()
    candidates = find_lib('nvvm', path)
    path = max(candidates) if candidates else None
    return _env_path_tuple(by, path)


def get_cuda_paths():
    """Returns a dictionary mapping component names to a 2-tuple
    of (source_variable, info).

    The returned dictionary will have the following keys and infos:
    - "nvvm": file_path
    - "libdevice": List[Tuple[arch, file_path]]
    - "cudalib_dir": directory_path

    Note: The result of the function is cached.
    """
    # Check cache
    if hasattr(get_cuda_paths, '_cached_result'):
        return get_cuda_paths._cached_result
    else:
        # Not in cache
        d = {
            'nvvm': _get_nvvm_path(),
            'libdevice': _get_libdevice_paths(),
            'cudalib_dir': _get_cudalib_dir(),
            'static_cudalib_dir': _get_static_cudalib_dir(),
        }
        # Cache result
        get_cuda_paths._cached_result = d
        return d


def get_debian_pkg_libdevice():
    """
    Return the Debian NVIDIA Maintainers-packaged libdevice location, if it
    exists.
    """
    pkg_libdevice_location = '/usr/lib/nvidia-cuda-toolkit/libdevice'
    if not os.path.exists(pkg_libdevice_location):
        return None
    return pkg_libdevice_location
