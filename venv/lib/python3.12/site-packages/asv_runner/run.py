import json
import math
import pickle

from ._aux import set_cpu_affinity_from_params
from .benchmarks.mark import SkipNotImplemented
from .discovery import get_benchmark_from_name


def _run(args):
    """
    Runs a specified benchmark and writes the result to a file.

    #### Parameters
    **args** (`tuple`)
    : A tuple containing benchmark directory, benchmark id,
    parameters string, profile path, and result file path.

    #### Notes
    This function first loads the extra parameters and sets the
    CPU affinity based on them. It then creates a benchmark
    from the `benchmark_id`. If the benchmark has a setup
    cache key, it loads the cache from a file and inserts it
    into the benchmark parameters.

    Then, the function runs the setup for the benchmark. If
    the setup indicates that the benchmark should be skipped,
    it sets the result as `math.nan`. Otherwise, it runs the
    benchmark and profiles it if a `profile_path` is provided.
    After running the benchmark, it performs the teardown for
    the benchmark and writes the result to the `result_file`.

    The `args` tuple contains:
    - **benchmark_dir** (`str`)
    : The directory where the benchmarks are located.
    - **benchmark_id** (`str`)
    : The id of the benchmark to run.
    - **params_str** (`str`)
    : A string containing JSON-encoded extra parameters.
    - **profile_path** (`str`)
    : The path for profile data. "None" implies no profiling.
    - **result_file** (`str`)
    : The path to the file where the result should be written.
    """

    (benchmark_dir, benchmark_id, params_str, profile_path, result_file) = args

    extra_params = json.loads(params_str)
    set_cpu_affinity_from_params(extra_params)
    extra_params.pop("cpu_affinity", None)

    if profile_path == "None":
        profile_path = None

    benchmark = get_benchmark_from_name(
        benchmark_dir, benchmark_id, extra_params=extra_params
    )

    if benchmark.setup_cache_key is not None:
        with open("cache.pickle", "rb") as fd:
            cache = pickle.load(fd)
        if cache is not None:
            benchmark.insert_param(cache)

    skip = benchmark.do_setup()

    try:
        if skip:
            result = math.nan
        else:
            try:
                result = benchmark.do_run()
                if profile_path is not None:
                    benchmark.do_profile(profile_path)
            except SkipNotImplemented:
                # Still runs setup() though
                result = math.nan
    finally:
        benchmark.do_teardown()

    with open(result_file, "w") as fp:
        json.dump(result, fp)
