import re

from ._base import Benchmark, _get_first_attr


class TrackBenchmark(Benchmark):
    """
    Represents a single benchmark for tracking an arbitrary value.

    The TrackBenchmark class provides a benchmark type for tracking any arbitrary
    value that your code produces. This can be useful when you need to track a value
    that isn't related to time or memory usage.

    #### Attributes
    **name_regex** (`re.Pattern`)
    : The regular expression used to match the names of functions that should be
    considered as track benchmarks.

    **type** (`str`)
    : The type of the benchmark. The default type is "track".

    **unit** (`str`)
    : The unit of the value that's being tracked. By default, this is "unit".

    #### Methods
    **run(*param)**
    : Runs the benchmark function and returns its result.

    """

    name_regex = re.compile("^(Track[A-Z_].+)|(track_.+)$")

    def __init__(self, name, func, attr_sources):
        """
        Initializes a new instance of the TrackBenchmark class.

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
        self.type = _get_first_attr(attr_sources, "type", "track")
        self.unit = _get_first_attr(attr_sources, "unit", "unit")

    def run(self, *param):
        """
        Runs the benchmark function and returns its result.

        #### Parameters
        **param** (`tuple`)
        : The parameters to pass to the benchmark function.

        #### Returns
        **result**
        : The result of the benchmark function.
        """
        return self.func(*param)


export_as_benchmark = [TrackBenchmark]
