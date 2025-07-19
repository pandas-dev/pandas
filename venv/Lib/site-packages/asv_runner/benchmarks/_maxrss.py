import os
import sys

ON_PYPY = hasattr(sys, "pypy_version_info")

if sys.platform.startswith("win"):
    import ctypes.wintypes

    SIZE_T = ctypes.c_size_t

    class PROCESS_MEMORY_COUNTERS(ctypes.Structure):
        """
        The PROCESS_MEMORY_COUNTERS structure is used by the
        GetProcessMemoryInfo function to store performance information. It's
        used here to retrieve the peak working set size, which is the maximum
        amount of memory in the working set of the process at any point in time.
        """

        _fields_ = [
            ("cb", ctypes.wintypes.DWORD),
            ("PageFaultCount", ctypes.wintypes.DWORD),
            ("PeakWorkingSetSize", SIZE_T),
            ("WorkingSetSize", SIZE_T),
            ("QuotaPeakPagedPoolUsage", SIZE_T),
            ("QuotaPagedPoolUsage", SIZE_T),
            ("QuotaPeakNonPagedPoolUsage", SIZE_T),
            ("QuotaNonPagedPoolUsage", SIZE_T),
            ("PagefileUsage", SIZE_T),
            ("PeakPagefileUsage", SIZE_T),
        ]

    GetCurrentProcess = ctypes.windll.kernel32.GetCurrentProcess
    GetCurrentProcess.argtypes = []
    GetCurrentProcess.restype = ctypes.wintypes.HANDLE

    GetProcessMemoryInfo = ctypes.windll.psapi.GetProcessMemoryInfo
    GetProcessMemoryInfo.argtypes = (
        ctypes.wintypes.HANDLE,
        ctypes.POINTER(PROCESS_MEMORY_COUNTERS),
        ctypes.wintypes.DWORD,
    )
    GetProcessMemoryInfo.restype = ctypes.wintypes.BOOL

    def get_maxrss():
        """
        Returns the peak working set size for the current process. On Windows,
        the peak working set size is the maximum amount of physical memory
        used by the process.

        #### Returns
        **peak_working_set_size** (`int`)
        : The peak working set size for the current process.
        """
        proc_hnd = GetCurrentProcess()
        counters = PROCESS_MEMORY_COUNTERS()
        info = GetProcessMemoryInfo(
            proc_hnd, ctypes.byref(counters), ctypes.sizeof(counters)
        )
        if info == 0:
            raise ctypes.WinError()
        return counters.PeakWorkingSetSize

    # Determine correct DWORD_PTR type for current Python version (32 or 64 bit)
    if ctypes.sizeof(ctypes.c_void_p) == ctypes.sizeof(ctypes.c_uint64):
        DWORD_PTR = ctypes.c_uint64
    elif ctypes.sizeof(ctypes.c_void_p) == ctypes.sizeof(ctypes.c_uint32):
        DWORD_PTR = ctypes.c_uint32

    SetProcessAffinityMask = ctypes.windll.kernel32.SetProcessAffinityMask
    SetProcessAffinityMask.argtypes = [ctypes.wintypes.HANDLE, DWORD_PTR]
    SetProcessAffinityMask.restype = bool

    GetCurrentProcess = ctypes.windll.kernel32.GetCurrentProcess
    GetCurrentProcess.argtypes = []
    GetCurrentProcess.restype = ctypes.wintypes.HANDLE

    def set_cpu_affinity(affinity_list):
        """
        Set CPU affinity to CPUs listed (numbered 0...n-1). CPU affinity
        is about binding and unbinding a process to a physical CPU or a
        range of CPUs, so that the process in question uses only a subset
        of the available CPUs.

        #### Parameters
        **affinity_list** (`list`)
        : A list of CPU cores to which the current process will be bound.
        """
        mask = 0
        for num in affinity_list:
            mask |= 2**num

        # Pseudohandle, doesn't need to be closed
        handle = GetCurrentProcess()
        ok = SetProcessAffinityMask(handle, mask)
        if not ok:
            raise RuntimeError("SetProcessAffinityMask failed")

else:
    try:
        import resource

        # POSIX
        if sys.platform == "darwin":

            def get_maxrss():
                """
                Returns the peak resident set size for the current process. On macOS,
                the peak resident set size is the maximum amount of memory occupied by
                the process's resident set at any point in time.

                #### Returns
                **peak_resident_set_size** (`int`)
                : The peak resident set size for the current process.
                """
                # OSX getrusage returns maxrss in bytes
                # https://developer.apple.com/library/mac/documentation/Darwin/Reference/ManPages/man2/getrusage.2.html
                return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss

        else:

            def get_maxrss():
                """
                Returns the peak resident set size for the current process. On Linux,
                the peak resident set size is the maximum amount of memory occupied by
                the process's resident set at any point in time.

                #### Returns
                **peak_resident_set_size** (`int`)
                : The peak resident set size for the current process.
                """
                # Linux, *BSD return maxrss in kilobytes
                return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss * 1024

    except ImportError:
        pass

    def set_cpu_affinity(affinity_list):
        """
        Set CPU affinity to CPUs listed (numbered 0...n-1). CPU affinity
        is about binding and unbinding a process to a physical CPU or a
        range of CPUs, so that the process in question uses only a subset
        of the available CPUs.

        #### Parameters
        **affinity_list** (`list`)
        : A list of CPU cores to which the current process will be bound.
        """
        if hasattr(os, "sched_setaffinity"):
            os.sched_setaffinity(0, affinity_list)
        else:
            import psutil

            p = psutil.Process()
            if hasattr(p, "cpu_affinity"):
                p.cpu_affinity(affinity_list)
