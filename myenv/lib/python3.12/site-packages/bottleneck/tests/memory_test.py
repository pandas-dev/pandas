import numpy as np
import sys
import bottleneck as bn
import pytest


@pytest.mark.skipif(
    sys.platform.startswith("win"), reason="resource module not available on windows"
)
def test_memory_leak():
    import resource

    arr = np.arange(1).reshape((1, 1))

    starting = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss

    for i in range(1000):
        for axis in [None, 0, 1]:
            bn.nansum(arr, axis=axis)
            bn.nanargmax(arr, axis=axis)
            bn.nanargmin(arr, axis=axis)
            bn.nanmedian(arr, axis=axis)
            bn.nansum(arr, axis=axis)
            bn.nanmean(arr, axis=axis)
            bn.nanmin(arr, axis=axis)
            bn.nanmax(arr, axis=axis)
            bn.nanvar(arr, axis=axis)

    ending = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss

    diff = ending - starting
    diff_bytes = diff * resource.getpagesize()
    print(diff_bytes)
    # For 1.3.0 release, this had value of ~100kB
    assert diff_bytes == 0
