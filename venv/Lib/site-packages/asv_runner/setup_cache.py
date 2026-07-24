import json
import pickle

from ._aux import set_cpu_affinity_from_params
from .discovery import get_benchmark_from_name


def _setup_cache(args):
    """
    Sets up a cache for a benchmark and pickles it into a file.

    #### Parameters
    **args** (`tuple`)
    : A tuple containing the benchmark directory, benchmark ID, and parameter
    string.

    - `benchmark_dir` (`str`): The directory where the benchmarks are located.
    - `benchmark_id` (`str`): The ID of the specific benchmark to be set up.
    - `params_str` (`str`): A JSON string containing extra parameters for the
      benchmark.

    #### Notes
    This function sets up a cache for a specific benchmark and saves it into a
    pickle file.

    First, it reads the extra parameters from the parameter string and sets the
    CPU affinity accordingly. Then, it retrieves the benchmark from the provided
    benchmark ID, and sets up its cache. The cache is then saved into a pickle
    file called 'cache.pickle'.

    The function is used to prepare a cache for a benchmark in a separate
    process, before running the actual benchmark.
    """
    (benchmark_dir, benchmark_id, params_str) = args

    extra_params = json.loads(params_str)

    set_cpu_affinity_from_params(extra_params)

    benchmark = get_benchmark_from_name(benchmark_dir, benchmark_id)
    cache = benchmark.do_setup_cache()
    with open("cache.pickle", "wb") as fd:
        pickle.dump(cache, fd)
