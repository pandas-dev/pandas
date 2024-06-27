import re

from ._base import Benchmark
from ._maxrss import get_maxrss


class PeakMemBenchmark(Benchmark):
    """
    Represents a single benchmark for tracking the peak memory consumption
    of the whole program.

    The PeakMemBenchmark class provides a benchmark type for tracking the peak
    memory consumption of the program while the benchmark function is running.

    #### Attributes
    **name_regex** (`re.Pattern`)
    : The regular expression used to match the names of functions that should be
    considered as peak memory benchmarks.

    **type** (`str`)
    : The type of the benchmark. The default type is "peakmemory".

    **unit** (`str`)
    : The unit of the value that's being tracked. By default, this is "bytes".

    #### Methods
    **run(*param)**
    : Runs the benchmark function and returns its result.

    """

    name_regex = re.compile("^(PeakMem[A-Z_].+)|(peakmem_.+)$")

    def __init__(self, name, func, attr_sources):
        """
        Initializes a new instance of the PeakMemBenchmark class.

        #### Parameters
        **name** (`str`)
        : The name of the benchmark.

        **func** (`callable`)
        : The function to benchmark.

        **attr_sources** (`list`)
        : A list of objects to search for attributes that might be used by the
        benchmark.
        """
        Benchmark.__init__(self, name, func, attr_sources)
        self.type = "peakmemory"
        self.unit = "bytes"

    def run(self, *param):
        """
        Runs the benchmark function and measures its peak memory consumption.

        #### Parameters
        **param** (`tuple`)
        : The parameters to pass to the benchmark function.

        #### Returns
        **result** (`int`)
        : The peak memory consumption in bytes of the program while the
        benchmark function was running.
        """
        self.func(*param)
        return get_maxrss()


export_as_benchmark = [PeakMemBenchmark]
