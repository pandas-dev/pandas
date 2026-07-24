import numpy as np
import bottleneck as bn
from .autotimeit import autotimeit

__all__ = ["bench_detailed"]


def bench_detailed(function="nansum", fraction_nan=0.0):
    """
    Benchmark a single function in detail or, optionally, all functions.

    Parameters
    ----------
    function : str, optional
        Name of function, as a string, to benchmark. Default ('nansum') is
        to benchmark bn.nansum. If `function` is 'all' then detailed
        benchmarks are run on all bottleneck functions.
    fraction_nan : float, optional
        Fraction of array elements that should, on average, be NaN. The
        default (0.0) is not to set any elements to NaN.

    Returns
    -------
    A benchmark report is printed to stdout.

    """

    if function == "all":
        # benchmark all bottleneck functions
        funcs = bn.get_functions("all", as_string=True)
        funcs.sort()
        for func in funcs:
            bench_detailed(func, fraction_nan)

    if fraction_nan < 0 or fraction_nan > 1:
        raise ValueError("`fraction_nan` must be between 0 and 1, inclusive")

    tab = "    "

    # Header
    print("%s benchmark" % function)
    print("%sBottleneck %s; Numpy %s" % (tab, bn.__version__, np.__version__))
    print("%sSpeed is NumPy time divided by Bottleneck time" % tab)
    if fraction_nan == 0:
        print("%sNone of the array elements are NaN" % tab)
    else:
        print(
            "%s%.1f%% of the array elements are NaN (on average)"
            % (tab, fraction_nan * 100)
        )
    print("")

    print("   Speed  Call                          Array")
    suite = benchsuite(function, fraction_nan)
    for test in suite:
        name = test["name"]
        speed = timer(test["statements"], test["setup"], test["repeat"])
        print("%8.1f  %s   %s" % (speed, name[0].ljust(27), name[1]))


def timer(statements, setup, repeat):
    if len(statements) != 2:
        raise ValueError("Two statements needed.")
    with np.errstate(invalid="ignore"):
        t0 = autotimeit(statements[0], setup, repeat=repeat)
        t1 = autotimeit(statements[1], setup, repeat=repeat)
    speed = t1 / t0
    return speed


def benchsuite(function, fraction_nan):

    # setup is called before each run of each function
    setup = """
        from bottleneck import %s as bn_fn
        try: from numpy import %s as sl_fn
        except ImportError: from bottleneck.slow import %s as sl_fn

        # avoid all-nan slice warnings from np.median and np.nanmedian
        if "%s" == "median": from bottleneck.slow import median as sl_fn
        if "%s" == "nanmedian": from bottleneck.slow import nanmedian as sl_fn

        from numpy import array, nan
        from numpy.random import RandomState
        rand = RandomState(123).rand

        a = %s
        if %s != 0: a[a < %s] = nan
    """
    setup = "\n".join([s.strip() for s in setup.split("\n")])

    # what kind of function signature do we need to use?
    if function in bn.get_functions("reduce", as_string=True):
        index = 0
    elif function in ["rankdata", "nanrankdata"]:
        index = 0
    elif function in bn.get_functions("move", as_string=True):
        index = 1
    elif function in ["partition", "argpartition", "push"]:
        index = 2
    elif function == "replace":
        index = 3
    else:
        raise ValueError("`function` (%s) not recognized" % function)

    # create benchmark suite
    instructions = get_instructions()
    f = function
    suite = []
    for instruction in instructions:
        signature = instruction[index + 1]
        if signature is None:
            continue
        array = instruction[0]
        repeat = instruction[-1]
        run = {}
        run["name"] = [f + signature, array]
        run["statements"] = ["bn_fn" + signature, "sl_fn" + signature]
        run["setup"] = setup % (f, f, f, f, f, array, fraction_nan, fraction_nan)
        run["repeat"] = repeat
        suite.append(run)

    return suite


def get_instructions():

    instructions = [
        # 1d input array
        (
            "rand(1)",
            "(a)",  # reduce + (nan)rankdata
            "(a, 1)",  # move
            "(a, 0)",  # (arg)partition
            "(a, np.nan, 0)",  # replace
            10,
        ),
        ("rand(10)", "(a)", "(a, 2)", "(a, 2)", "(a, np.nan, 0)", 10),
        ("rand(100)", "(a)", "(a, 20)", "(a, 20)", "(a, np.nan, 0)", 6),
        ("rand(1000)", "(a)", "(a, 200)", "(a, 200)", "(a, np.nan, 0)", 3),
        ("rand(1000000)", "(a)", "(a, 200)", "(a, 200)", "(a, np.nan, 0)", 2),
        # 2d input array
        ("rand(10, 10)", "(a)", "(a, 2)", "(a, 2)", "(a, np.nan, 0)", 6),
        ("rand(100, 100)", "(a)", "(a, 20)", "(a, 20)", "(a, np.nan, 0)", 3),
        ("rand(1000, 1000)", "(a)", "(a, 200)", "(a, 200)", "(a, np.nan, 0)", 2),
        ("rand(10, 10)", "(a, 1)", None, None, None, 6),
        ("rand(100, 100)", "(a, 1)", None, None, None, 3),
        ("rand(1000, 1000)", "(a, 1)", None, None, None, 2),
        ("rand(100000, 2)", "(a, 1)", "(a, 1)", "(a, 1)", None, 2),
        ("rand(10, 10)", "(a, 0)", None, None, None, 6),
        ("rand(100, 100)", "(a, 0)", "(a, 20, axis=0)", None, None, 3),
        ("rand(1000, 1000)", "(a, 0)", "(a, 200, axis=0)", None, None, 2),
        # 3d input array
        (
            "rand(100, 100, 100)",
            "(a, 0)",
            "(a, 20, axis=0)",
            "(a, 20, axis=0)",
            None,
            2,
        ),
        (
            "rand(100, 100, 100)",
            "(a, 1)",
            "(a, 20, axis=1)",
            "(a, 20, axis=1)",
            None,
            2,
        ),
        (
            "rand(100, 100, 100)",
            "(a, 2)",
            "(a, 20, axis=2)",
            "(a, 20, axis=2)",
            "(a, np.nan, 0)",
            2,
        ),
        # 0d input array
        ("array(1.0)", "(a)", None, None, "(a, 0, 2)", 10),
    ]

    return instructions
