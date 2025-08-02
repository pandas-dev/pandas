import copy
import re

from ._base import Benchmark
from ._exceptions import NotRequired

try:
    from pympler.asizeof import asizeof
except ImportError:
    raise NotRequired("MemBenchmarks not requested or pympler not found")


class MemBenchmark(Benchmark):
    """
    Represents a single benchmark for tracking the memory consumption of an object.

    The MemBenchmark class provides a benchmark type for tracking the memory
    consumption of the object returned by the benchmark function.

    #### Attributes
    **name_regex** (`re.Pattern`)
    : The regular expression used to match the names of functions that should be
    considered as memory benchmarks.

    **type** (`str`)
    : The type of the benchmark. The default type is "memory".

    **unit** (`str`)
    : The unit of the value that's being tracked. By default, this is "bytes".

    #### Methods
    **run(*param)**
    : Runs the benchmark function and returns the memory consumption of the object
    returned by the function.

    """

    name_regex = re.compile("^(Mem[A-Z_].+)|(mem_.+)$")

    def __init__(self, name, func, attr_sources):
        """
        Initializes a new instance of the MemBenchmark class.

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
        self.type = "memory"
        self.unit = "bytes"

    def run(self, *param):
        """
        Runs the benchmark function and measures the memory consumption of the object
        returned by the function.

        #### Parameters
        **param** (`tuple`)
        : The parameters to pass to the benchmark function.

        #### Returns
        **result** (`int`)
        : The memory consumption in bytes of the object returned by the
        benchmark function.
        """
        obj = self.func(*param)

        sizeof2 = asizeof([obj, obj])
        sizeofcopy = asizeof([obj, copy.copy(obj)])

        return sizeofcopy - sizeof2


export_as_benchmark = [MemBenchmark]
