import json
import locale
import multiprocessing
import os
import platform
import textwrap
import sys
from contextlib import redirect_stdout
from datetime import datetime
from io import StringIO
from subprocess import check_output, PIPE, CalledProcessError
import numpy as np
import llvmlite.binding as llvmbind
from llvmlite import __version__ as llvmlite_version
from numba import cuda as cu, __version__ as version_number
from numba.cuda import cudadrv
from numba.cuda.cudadrv.driver import driver as cudriver
from numba.cuda.cudadrv.runtime import runtime as curuntime
from numba.core import config

_psutil_import = False
try:
    import psutil
except ImportError:
    pass
else:
    _psutil_import = True

__all__ = ['get_sysinfo', 'display_sysinfo']

# Keys of a `sysinfo` dictionary

# Time info
_start, _start_utc, _runtime = 'Start', 'Start UTC', 'Runtime'
_numba_version = 'Numba Version'
# Hardware info
_machine = 'Machine'
_cpu_name, _cpu_count = 'CPU Name', 'CPU Count'
_cpus_allowed, _cpus_list = 'CPUs Allowed', 'List CPUs Allowed'
_cpu_features = 'CPU Features'
_cfs_quota, _cfs_period = 'CFS Quota', 'CFS Period',
_cfs_restrict = 'CFS Restriction'
_mem_total, _mem_available = 'Mem Total', 'Mem Available'
# OS info
_platform_name, _platform_release = 'Platform Name', 'Platform Release'
_os_name, _os_version = 'OS Name', 'OS Version'
_os_spec_version = 'OS Specific Version'
_libc_version = 'Libc Version'
# Python info
_python_comp = 'Python Compiler'
_python_impl = 'Python Implementation'
_python_version = 'Python Version'
_python_locale = 'Python Locale'
# LLVM info
_llvmlite_version = 'llvmlite Version'
_llvm_version = 'LLVM Version'
# CUDA info
_cu_target_impl = 'CUDA Target Impl'
_cu_dev_init = 'CUDA Device Init'
_cu_drv_ver = 'CUDA Driver Version'
_cu_rt_ver = 'CUDA Runtime Version'
_cu_nvidia_bindings = 'NVIDIA CUDA Bindings'
_cu_nvidia_bindings_used = 'NVIDIA CUDA Bindings In Use'
_cu_detect_out, _cu_lib_test = 'CUDA Detect Output', 'CUDA Lib Test'
_cu_mvc_available = 'NVIDIA CUDA Minor Version Compatibility Available'
_cu_mvc_needed = 'NVIDIA CUDA Minor Version Compatibility Needed'
_cu_mvc_in_use = 'NVIDIA CUDA Minor Version Compatibility In Use'
# NumPy info
_numpy_version = 'NumPy Version'
_numpy_supported_simd_features = 'NumPy Supported SIMD features'
_numpy_supported_simd_dispatch = 'NumPy Supported SIMD dispatch'
_numpy_supported_simd_baseline = 'NumPy Supported SIMD baseline'
_numpy_AVX512_SKX_detected = 'NumPy AVX512_SKX detected'
# SVML info
_svml_state, _svml_loaded = 'SVML State', 'SVML Lib Loaded'
_llvm_svml_patched = 'LLVM SVML Patched'
_svml_operational = 'SVML Operational'
# Threading layer info
_tbb_thread, _tbb_error = 'TBB Threading', 'TBB Threading Error'
_openmp_thread, _openmp_error = 'OpenMP Threading', 'OpenMP Threading Error'
_openmp_vendor = 'OpenMP vendor'
_wkq_thread, _wkq_error = 'Workqueue Threading', 'Workqueue Threading Error'
# Numba info
_numba_env_vars = 'Numba Env Vars'
# Conda info
_conda_build_ver, _conda_env_ver = 'Conda Build', 'Conda Env'
_conda_platform, _conda_python_ver = 'Conda Platform', 'Conda Python Version'
_conda_root_writable = 'Conda Root Writable'
# Packages info
_inst_pkg = 'Installed Packages'
# Psutil info
_psutil = 'Psutil Available'
# Errors and warnings
_errors = 'Errors'
_warnings = 'Warnings'

# Error and warning log
_error_log = []
_warning_log = []


def get_os_spec_info(os_name):
    # Linux man page for `/proc`:
    # http://man7.org/linux/man-pages/man5/proc.5.html

    # Windows documentation for `wmic OS`:
    # https://docs.microsoft.com/en-us/windows/win32/cimwin32prov/cim-operatingsystem

    # MacOS man page for `sysctl`:
    # https://www.unix.com/man-page/osx/3/sysctl/
    # MacOS man page for `vm_stat`:
    # https://www.unix.com/man-page/osx/1/vm_stat/

    class CmdBufferOut(tuple):
        buffer_output_flag = True

    class CmdReadFile(tuple):
        read_file_flag = True

    shell_params = {
        'Linux': {
            'cmd': (
                CmdReadFile(('/sys/fs/cgroup/cpuacct/cpu.cfs_quota_us',)),
                CmdReadFile(('/sys/fs/cgroup/cpuacct/cpu.cfs_period_us',)),
            ),
            'cmd_optional': (
                CmdReadFile(('/proc/meminfo',)),
                CmdReadFile(('/proc/self/status',)),
            ),
            'kwds': {
                # output string fragment -> result dict key
                'MemTotal:': _mem_total,
                'MemAvailable:': _mem_available,
                'Cpus_allowed:': _cpus_allowed,
                'Cpus_allowed_list:': _cpus_list,
                '/sys/fs/cgroup/cpuacct/cpu.cfs_quota_us': _cfs_quota,
                '/sys/fs/cgroup/cpuacct/cpu.cfs_period_us': _cfs_period,
            },
        },
        'Windows': {
            'cmd': (),
            'cmd_optional': (
                CmdBufferOut(('wmic', 'OS', 'get', 'TotalVirtualMemorySize')),
                CmdBufferOut(('wmic', 'OS', 'get', 'FreeVirtualMemory')),
            ),
            'kwds': {
                # output string fragment -> result dict key
                'TotalVirtualMemorySize': _mem_total,
                'FreeVirtualMemory': _mem_available,
            },
        },
        'Darwin': {
            'cmd': (),
            'cmd_optional': (
                ('sysctl', 'hw.memsize'),
                ('vm_stat'),
            ),
            'kwds': {
                # output string fragment -> result dict key
                'hw.memsize:': _mem_total,
                'free:': _mem_available,
            },
            'units': {
                _mem_total: 1,  # Size is given in bytes.
                _mem_available: 4096,  # Size is given in 4kB pages.
            },
        },
    }

    os_spec_info = {}
    params = shell_params.get(os_name, {})
    cmd_selected = params.get('cmd', ())

    if _psutil_import:
        vm = psutil.virtual_memory()
        os_spec_info.update({
            _mem_total: vm.total,
            _mem_available: vm.available,
        })
        p = psutil.Process()
        cpus_allowed = p.cpu_affinity() if hasattr(p, 'cpu_affinity') else []
        if cpus_allowed:
            os_spec_info[_cpus_allowed] = len(cpus_allowed)
            os_spec_info[_cpus_list] = ' '.join(str(n) for n in cpus_allowed)

    else:
        _warning_log.append(
            "Warning (psutil): psutil cannot be imported. "
            "For more accuracy, consider installing it.")
        # Fallback to internal heuristics
        cmd_selected += params.get('cmd_optional', ())

    # Assuming the shell cmd returns a unique (k, v) pair per line
    # or a unique (k, v) pair spread over several lines:
    # Gather output in a list of strings containing a keyword and some value.
    output = []
    for cmd in cmd_selected:
        if hasattr(cmd, 'read_file_flag'):
            # Open file within Python
            if os.path.exists(cmd[0]):
                try:
                    with open(cmd[0], 'r') as f:
                        out = f.readlines()
                        if out:
                            out[0] = ' '.join((cmd[0], out[0]))
                            output.extend(out)
                except OSError as e:
                    _error_log.append(f'Error (file read): {e}')
                    continue
            else:
                _warning_log.append('Warning (no file): {}'.format(cmd[0]))
                continue
        else:
            # Spawn a subprocess
            try:
                out = check_output(cmd, stderr=PIPE)
            except (OSError, CalledProcessError) as e:
                _error_log.append(f'Error (subprocess): {e}')
                continue
            if hasattr(cmd, 'buffer_output_flag'):
                out = b' '.join(line for line in out.splitlines()) + b'\n'
            output.extend(out.decode().splitlines())

    # Extract (k, output) pairs by searching for keywords in output
    kwds = params.get('kwds', {})
    for line in output:
        match = kwds.keys() & line.split()
        if match and len(match) == 1:
            k = kwds[match.pop()]
            os_spec_info[k] = line
        elif len(match) > 1:
            print(f'Ambiguous output: {line}')

    # Try to extract something meaningful from output string
    def format():
        # CFS restrictions
        split = os_spec_info.get(_cfs_quota, '').split()
        if split:
            os_spec_info[_cfs_quota] = float(split[-1])
        split = os_spec_info.get(_cfs_period, '').split()
        if split:
            os_spec_info[_cfs_period] = float(split[-1])
        if os_spec_info.get(_cfs_quota, -1) != -1:
            cfs_quota = os_spec_info.get(_cfs_quota, '')
            cfs_period = os_spec_info.get(_cfs_period, '')
            runtime_amount = cfs_quota / cfs_period
            os_spec_info[_cfs_restrict] = runtime_amount

    def format_optional():
        # Memory
        units = {_mem_total: 1024, _mem_available: 1024}
        units.update(params.get('units', {}))
        for k in (_mem_total, _mem_available):
            digits = ''.join(d for d in os_spec_info.get(k, '') if d.isdigit())
            os_spec_info[k] = int(digits or 0) * units[k]
        # Accessible CPUs
        split = os_spec_info.get(_cpus_allowed, '').split()
        if split:
            n = split[-1]
            n = n.split(',')[-1]
            os_spec_info[_cpus_allowed] = str(bin(int(n or 0, 16))).count('1')
        split = os_spec_info.get(_cpus_list, '').split()
        if split:
            os_spec_info[_cpus_list] = split[-1]

    try:
        format()
        if not _psutil_import:
            format_optional()
    except Exception as e:
        _error_log.append(f'Error (format shell output): {e}')

    # Call OS specific functions
    os_specific_funcs = {
        'Linux': {
            _libc_version: lambda: ' '.join(platform.libc_ver())
        },
        'Windows': {
            _os_spec_version: lambda: ' '.join(
                s for s in platform.win32_ver()),
        },
        'Darwin': {
            _os_spec_version: lambda: ''.join(
                i or ' ' for s in tuple(platform.mac_ver()) for i in s),
        },
    }
    key_func = os_specific_funcs.get(os_name, {})
    os_spec_info.update({k: f() for k, f in key_func.items()})
    return os_spec_info


def get_sysinfo():

    # Gather the information that shouldn't raise exceptions
    sys_info = {
        _start: datetime.now(),
        _start_utc: datetime.utcnow(),
        _machine: platform.machine(),
        _cpu_name: llvmbind.get_host_cpu_name(),
        _cpu_count: multiprocessing.cpu_count(),
        _platform_name: platform.platform(aliased=True),
        _platform_release: platform.release(),
        _os_name: platform.system(),
        _os_version: platform.version(),
        _python_comp: platform.python_compiler(),
        _python_impl: platform.python_implementation(),
        _python_version: platform.python_version(),
        _numba_env_vars: {k: v for (k, v) in os.environ.items()
                          if k.startswith('NUMBA_')},
        _numba_version: version_number,
        _llvm_version: '.'.join(str(i) for i in llvmbind.llvm_version_info),
        _llvmlite_version: llvmlite_version,
        _psutil: _psutil_import,
    }

    # CPU features
    try:
        feature_map = llvmbind.get_host_cpu_features()
    except RuntimeError as e:
        _error_log.append(f'Error (CPU features): {e}')
    else:
        features = sorted([key for key, value in feature_map.items() if value])
        sys_info[_cpu_features] = ' '.join(features)

    # Python locale
    # On MacOSX, getdefaultlocale can raise. Check again if Py > 3.7.5
    try:
        # If $LANG is unset, getdefaultlocale() can return (None, None), make
        # sure we can encode this as strings by casting explicitly.
        sys_info[_python_locale] = '.'.join([str(i) for i in
                                             locale.getdefaultlocale()])
    except Exception as e:
        _error_log.append(f'Error (locale): {e}')

    # CUDA information
    try:
        sys_info[_cu_target_impl] = cu.implementation
    except AttributeError:
        # On the offchance an out-of-tree target did not set the
        # implementation, we can try to continue
        pass

    try:
        cu.list_devices()[0]  # will a device initialise?
    except Exception as e:
        sys_info[_cu_dev_init] = False
        msg_not_found = "CUDA driver library cannot be found"
        msg_disabled_by_user = "CUDA is disabled"
        msg_end = " or no CUDA enabled devices are present."
        msg_generic_problem = "CUDA device initialisation problem."
        msg = getattr(e, 'msg', None)
        if msg is not None:
            if msg_not_found in msg:
                err_msg = msg_not_found + msg_end
            elif msg_disabled_by_user in msg:
                err_msg = msg_disabled_by_user + msg_end
            else:
                err_msg = msg_generic_problem + " Message:" + msg
        else:
            err_msg = msg_generic_problem + " " + str(e)
        # Best effort error report
        _warning_log.append("Warning (cuda): %s\nException class: %s" %
                            (err_msg, str(type(e))))
    else:
        try:
            sys_info[_cu_dev_init] = True

            output = StringIO()
            with redirect_stdout(output):
                cu.detect()
            sys_info[_cu_detect_out] = output.getvalue()
            output.close()

            cu_drv_ver = cudriver.get_version()
            cu_rt_ver = curuntime.get_version()
            sys_info[_cu_drv_ver] = '%s.%s' % cu_drv_ver
            sys_info[_cu_rt_ver] = '%s.%s' % cu_rt_ver

            output = StringIO()
            with redirect_stdout(output):
                cudadrv.libs.test()
            sys_info[_cu_lib_test] = output.getvalue()
            output.close()

            try:
                from cuda import cuda  # noqa: F401
                nvidia_bindings_available = True
            except ImportError:
                nvidia_bindings_available = False
            sys_info[_cu_nvidia_bindings] = nvidia_bindings_available

            nv_binding_used = bool(cudadrv.driver.USE_NV_BINDING)
            sys_info[_cu_nvidia_bindings_used] = nv_binding_used

            try:
                from ptxcompiler import compile_ptx  # noqa: F401
                from cubinlinker import CubinLinker  # noqa: F401
                sys_info[_cu_mvc_available] = True
            except ImportError:
                sys_info[_cu_mvc_available] = False

            sys_info[_cu_mvc_needed] = cu_rt_ver > cu_drv_ver
            sys_info[_cu_mvc_in_use] = bool(
                config.CUDA_ENABLE_MINOR_VERSION_COMPATIBILITY)
        except Exception as e:
            _warning_log.append(
                "Warning (cuda): Probing CUDA failed "
                "(device and driver present, runtime problem?)\n"
                f"(cuda) {type(e)}: {e}")

    # NumPy information
    sys_info[_numpy_version] = np.version.full_version
    try:
        # NOTE: These consts were added in NumPy 1.20
        from numpy.core._multiarray_umath import (__cpu_features__,
                                                  __cpu_dispatch__,
                                                  __cpu_baseline__,)
    except ImportError:
        sys_info[_numpy_AVX512_SKX_detected] = False
    else:
        feat_filtered = [k for k, v in __cpu_features__.items() if v]
        sys_info[_numpy_supported_simd_features] = feat_filtered
        sys_info[_numpy_supported_simd_dispatch] = __cpu_dispatch__
        sys_info[_numpy_supported_simd_baseline] = __cpu_baseline__
        sys_info[_numpy_AVX512_SKX_detected] = \
            __cpu_features__.get("AVX512_SKX", False)

    # SVML information
    # Replicate some SVML detection logic from numba.__init__ here.
    # If SVML load fails in numba.__init__ the splitting of the logic
    # here will help diagnosing the underlying issue.
    svml_lib_loaded = True
    try:
        if sys.platform.startswith('linux'):
            llvmbind.load_library_permanently("libsvml.so")
        elif sys.platform.startswith('darwin'):
            llvmbind.load_library_permanently("libsvml.dylib")
        elif sys.platform.startswith('win'):
            llvmbind.load_library_permanently("svml_dispmd")
        else:
            svml_lib_loaded = False
    except Exception:
        svml_lib_loaded = False
    func = getattr(llvmbind.targets, "has_svml", None)
    sys_info[_llvm_svml_patched] = func() if func else False
    sys_info[_svml_state] = config.USING_SVML
    sys_info[_svml_loaded] = svml_lib_loaded
    sys_info[_svml_operational] = all((
        sys_info[_svml_state],
        sys_info[_svml_loaded],
        sys_info[_llvm_svml_patched],
    ))

    # Check which threading backends are available.
    def parse_error(e, backend):
        # parses a linux based error message, this is to provide feedback
        # and hide user paths etc
        try:
            path, problem, symbol = [x.strip() for x in e.msg.split(':')]
            extn_dso = os.path.split(path)[1]
            if backend in extn_dso:
                return "%s: %s" % (problem, symbol)
        except Exception:
            pass
        return "Unknown import problem."

    try:
        # check import is ok, this means the DSO linkage is working
        from numba.np.ufunc import tbbpool  # NOQA
        # check that the version is compatible, this is a check performed at
        # runtime (well, compile time), it will also ImportError if there's
        # a problem.
        from numba.np.ufunc.parallel import _check_tbb_version_compatible
        _check_tbb_version_compatible()
        sys_info[_tbb_thread] = True
    except ImportError as e:
        # might be a missing symbol due to e.g. tbb libraries missing
        sys_info[_tbb_thread] = False
        sys_info[_tbb_error] = parse_error(e, 'tbbpool')

    try:
        from numba.np.ufunc import omppool
        sys_info[_openmp_thread] = True
        sys_info[_openmp_vendor] = omppool.openmp_vendor
    except ImportError as e:
        sys_info[_openmp_thread] = False
        sys_info[_openmp_error] = parse_error(e, 'omppool')

    try:
        from numba.np.ufunc import workqueue  # NOQA
        sys_info[_wkq_thread] = True
    except ImportError as e:
        sys_info[_wkq_thread] = True
        sys_info[_wkq_error] = parse_error(e, 'workqueue')

    # Look for conda and installed packages information
    cmd = ('conda', 'info', '--json')
    try:
        conda_out = check_output(cmd)
    except Exception as e:
        _warning_log.append(f'Warning: Conda not available.\n Error was {e}\n')
        # Conda is not available, try pip list to list installed packages
        cmd = (sys.executable, '-m', 'pip', 'list')
        try:
            reqs = check_output(cmd)
        except Exception as e:
            _error_log.append(f'Error (pip): {e}')
        else:
            sys_info[_inst_pkg] = reqs.decode().splitlines()

    else:
        jsond = json.loads(conda_out.decode())
        keys = {
            'conda_build_version': _conda_build_ver,
            'conda_env_version': _conda_env_ver,
            'platform': _conda_platform,
            'python_version': _conda_python_ver,
            'root_writable': _conda_root_writable,
        }
        for conda_k, sysinfo_k in keys.items():
            sys_info[sysinfo_k] = jsond.get(conda_k, 'N/A')

        # Get info about packages in current environment
        cmd = ('conda', 'list')
        try:
            conda_out = check_output(cmd)
        except CalledProcessError as e:
            _error_log.append(f'Error (conda): {e}')
        else:
            data = conda_out.decode().splitlines()
            sys_info[_inst_pkg] = [l for l in data if not l.startswith('#')]

    sys_info.update(get_os_spec_info(sys_info[_os_name]))
    sys_info[_errors] = _error_log
    sys_info[_warnings] = _warning_log
    sys_info[_runtime] = (datetime.now() - sys_info[_start]).total_seconds()
    return sys_info


def display_sysinfo(info=None, sep_pos=45):
    class DisplayMap(dict):
        display_map_flag = True

    class DisplaySeq(tuple):
        display_seq_flag = True

    class DisplaySeqMaps(tuple):
        display_seqmaps_flag = True

    if info is None:
        info = get_sysinfo()

    fmt = f'%-{sep_pos}s : %-s'
    MB = 1024**2
    template = (
        ("-" * 80,),
        ("__Time Stamp__",),
        ("Report started (local time)", info.get(_start, '?')),
        ("UTC start time", info.get(_start_utc, '?')),
        ("Running time (s)", info.get(_runtime, '?')),
        ("",),
        ("__Hardware Information__",),
        ("Machine", info.get(_machine, '?')),
        ("CPU Name", info.get(_cpu_name, '?')),
        ("CPU Count", info.get(_cpu_count, '?')),
        ("Number of accessible CPUs", info.get(_cpus_allowed, '?')),
        ("List of accessible CPUs cores", info.get(_cpus_list, '?')),
        ("CFS Restrictions (CPUs worth of runtime)",
            info.get(_cfs_restrict, 'None')),
        ("",),
        ("CPU Features", '\n'.join(
            ' ' * (sep_pos + 3) + l if i else l
            for i, l in enumerate(
                textwrap.wrap(
                    info.get(_cpu_features, '?'),
                    width=79 - sep_pos
                )
            )
        )),
        ("",),
        ("Memory Total (MB)", info.get(_mem_total, 0) // MB or '?'),
        ("Memory Available (MB)"
            if info.get(_os_name, '') != 'Darwin' or info.get(_psutil, False)
            else "Free Memory (MB)", info.get(_mem_available, 0) // MB or '?'),
        ("",),
        ("__OS Information__",),
        ("Platform Name", info.get(_platform_name, '?')),
        ("Platform Release", info.get(_platform_release, '?')),
        ("OS Name", info.get(_os_name, '?')),
        ("OS Version", info.get(_os_version, '?')),
        ("OS Specific Version", info.get(_os_spec_version, '?')),
        ("Libc Version", info.get(_libc_version, '?')),
        ("",),
        ("__Python Information__",),
        DisplayMap({k: v for k, v in info.items() if k.startswith('Python')}),
        ("",),
        ("__Numba Toolchain Versions__",),
        ("Numba Version", info.get(_numba_version, '?')),
        ("llvmlite Version", info.get(_llvmlite_version, '?')),
        ("",),
        ("__LLVM Information__",),
        ("LLVM Version", info.get(_llvm_version, '?')),
        ("",),
        ("__CUDA Information__",),
        ("CUDA Target Implementation", info.get(_cu_target_impl, '?')),
        ("CUDA Device Initialized", info.get(_cu_dev_init, '?')),
        ("CUDA Driver Version", info.get(_cu_drv_ver, '?')),
        ("CUDA Runtime Version", info.get(_cu_rt_ver, '?')),
        ("CUDA NVIDIA Bindings Available", info.get(_cu_nvidia_bindings, '?')),
        ("CUDA NVIDIA Bindings In Use",
         info.get(_cu_nvidia_bindings_used, '?')),
        ("CUDA Minor Version Compatibility Available",
         info.get(_cu_mvc_available, '?')),
        ("CUDA Minor Version Compatibility Needed",
         info.get(_cu_mvc_needed, '?')),
        ("CUDA Minor Version Compatibility In Use",
         info.get(_cu_mvc_in_use, '?')),
        ("CUDA Detect Output:",),
        (info.get(_cu_detect_out, "None"),),
        ("CUDA Libraries Test Output:",),
        (info.get(_cu_lib_test, "None"),),
        ("",),
        ("__NumPy Information__",),
        ("NumPy Version", info.get(_numpy_version, '?')),
        ("NumPy Supported SIMD features",
         DisplaySeq(info.get(_numpy_supported_simd_features, [])
                    or ('None found.',))),
        ("NumPy Supported SIMD dispatch",
         DisplaySeq(info.get(_numpy_supported_simd_dispatch, [])
                    or ('None found.',))),
        ("NumPy Supported SIMD baseline",
         DisplaySeq(info.get(_numpy_supported_simd_baseline, [])
                    or ('None found.',))),
        ("NumPy AVX512_SKX support detected",
         info.get(_numpy_AVX512_SKX_detected, '?')),
        ("",),
        ("__SVML Information__",),
        ("SVML State, config.USING_SVML", info.get(_svml_state, '?')),
        ("SVML Library Loaded", info.get(_svml_loaded, '?')),
        ("llvmlite Using SVML Patched LLVM", info.get(_llvm_svml_patched, '?')),
        ("SVML Operational", info.get(_svml_operational, '?')),
        ("",),
        ("__Threading Layer Information__",),
        ("TBB Threading Layer Available", info.get(_tbb_thread, '?')),
        ("+-->TBB imported successfully." if info.get(_tbb_thread, '?')
            else f"+--> Disabled due to {info.get(_tbb_error, '?')}",),
        ("OpenMP Threading Layer Available", info.get(_openmp_thread, '?')),
        (f"+-->Vendor: {info.get(_openmp_vendor, '?')}"
            if info.get(_openmp_thread, False)
            else f"+--> Disabled due to {info.get(_openmp_error, '?')}",),
        ("Workqueue Threading Layer Available", info.get(_wkq_thread, '?')),
        ("+-->Workqueue imported successfully." if info.get(_wkq_thread, False)
            else f"+--> Disabled due to {info.get(_wkq_error, '?')}",),
        ("",),
        ("__Numba Environment Variable Information__",),
        (DisplayMap(info.get(_numba_env_vars, {})) or ('None found.',)),
        ("",),
        ("__Conda Information__",),
        (DisplayMap({k: v for k, v in info.items()
                     if k.startswith('Conda')}) or ("Conda not available.",)),
        ("",),
        ("__Installed Packages__",),
        DisplaySeq(info.get(_inst_pkg, ("Couldn't retrieve packages info.",))),
        ("",),
        ("__Error log__" if info.get(_errors, [])
            else "No errors reported.",),
        DisplaySeq(info.get(_errors, [])),
        ("",),
        ("__Warning log__" if info.get(_warnings, [])
            else "No warnings reported.",),
        DisplaySeq(info.get(_warnings, [])),
        ("-" * 80,),
        ("If requested, please copy and paste the information between\n"
         "the dashed (----) lines, or from a given specific section as\n"
         "appropriate.\n\n"
         "=============================================================\n"
         "IMPORTANT: Please ensure that you are happy with sharing the\n"
         "contents of the information present, any information that you\n"
         "wish to keep private you should remove before sharing.\n"
         "=============================================================\n",),
    )
    for t in template:
        if hasattr(t, 'display_seq_flag'):
            print(*t, sep='\n')
        elif hasattr(t, 'display_map_flag'):
            print(*tuple(fmt % (k, v) for (k, v) in t.items()), sep='\n')
        elif hasattr(t, 'display_seqmaps_flag'):
            for d in t:
                print(*tuple(fmt % ('\t' + k, v) for (k, v) in d.items()),
                      sep='\n', end='\n')
        elif len(t) == 2:
            print(fmt % t)
        else:
            print(*t)


if __name__ == '__main__':
    display_sysinfo()
