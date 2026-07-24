import numpy as np
import bottleneck as bn
from .autotimeit import autotimeit

__all__ = ["bench"]


def bench(
    shapes=[(100,), (1000, 1000), (1000, 1000), (1000, 1000), (1000, 1000)],
    axes=[0, 0, 0, 1, 1],
    nans=[False, False, True, False, True],
    dtype="float64",
    order="C",
    functions=None,
):
    """
    Bottleneck benchmark.

    Parameters
    ----------
    shapes : list, optional
        A list of tuple shapes of input arrays to use in the benchmark.
    axes : list, optional
        List of axes along which to perform the calculations that are being
        benchmarked.
    nans : list, optional
        A list of the bools (True or False), one for each tuple in the
        `shapes` list, that tells whether the input arrays should be randomly
        filled with one-fifth NaNs.
    dtype : str, optional
        Data type string such as 'float64', which is the default.
    order : {'C', 'F'}, optional
        Whether to store multidimensional data in C- or Fortran-contiguous
        (row- or column-wise) order in memory.
    functions : {list, None}, optional
        A list of strings specifying which functions to include in the
        benchmark. By default (None) all functions are included in the
        benchmark.

    Returns
    -------
    A benchmark report is printed to stdout.

    """

    if len(shapes) != len(nans):
        raise ValueError("`shapes` and `nans` must have the same length")
    if len(shapes) != len(axes):
        raise ValueError("`shapes` and `axes` must have the same length")

    # Header
    print("Bottleneck performance benchmark")
    print("    Bottleneck %s; Numpy %s" % (bn.__version__, np.__version__))
    print("    Speed is NumPy time divided by Bottleneck time")
    print("    NaN means approx one-fifth NaNs; %s used" % str(dtype))

    print("")
    header = [" " * 11]
    for nan in nans:
        if nan:
            header.append("NaN".center(11))
        else:
            header.append("no NaN".center(11))
    print("".join(header))
    header = ["".join(str(shape).split(" ")).center(11) for shape in shapes]
    header = [" " * 12] + header
    print("".join(header))
    header = ["".join(("axis=" + str(axis)).split(" ")).center(11) for axis in axes]
    header = [" " * 12] + header
    print("".join(header))

    suite = benchsuite(shapes, dtype, nans, axes, order, functions)
    for test in suite:
        name = test["name"].ljust(12)
        fmt = name + "%7.1f" + "%11.1f" * (len(shapes) - 1)
        speed = timer(test["statements"], test["setups"])
        print(fmt % tuple(speed))


def timer(statements, setups):
    speed = []
    if len(statements) != 2:
        raise ValueError("Two statements needed.")
    for setup in setups:
        with np.errstate(invalid="ignore"):
            t0 = autotimeit(statements[0], setup)
            t1 = autotimeit(statements[1], setup)
        speed.append(t1 / t0)
    return speed


def getarray(shape, dtype, nans, order):
    a = np.arange(np.prod(shape), dtype=dtype)
    if nans and issubclass(a.dtype.type, np.inexact):
        a[::5] = np.nan
    rs = np.random.RandomState(shape)
    rs.shuffle(a)
    return np.array(a.reshape(*shape), order=order)


def benchsuite(shapes, dtype, nans, axes, order, functions):

    suite = []

    def getsetups(setup, shapes, nans, axes, dtype, order):
        template = """
        from bottleneck.benchmark.bench import getarray
        a = getarray(%s, '%s', %s, '%s')
        axis=%s
        %s"""
        setups = []
        for shape, axis, nan in zip(shapes, axes, nans):
            s = template % (
                str(shape),
                str(dtype),
                str(nan),
                str(order),
                str(axis),
                setup,
            )
            s = "\n".join([line.strip() for line in s.split("\n")])
            setups.append(s)
        return setups

    # non-moving window functions
    funcs = bn.get_functions("reduce", as_string=True)
    funcs += ["rankdata", "nanrankdata"]
    for func in funcs:
        if functions is not None and func not in functions:
            continue
        run = {}
        run["name"] = func
        run["statements"] = ["bn_func(a, axis)", "sl_func(a, axis)"]
        setup = """
            from bottleneck import %s as bn_func
            try: from numpy import %s as sl_func
            except ImportError: from bottleneck.slow import %s as sl_func
            if "%s" == "median": from bottleneck.slow import median as sl_func
        """ % (
            func,
            func,
            func,
            func,
        )
        run["setups"] = getsetups(setup, shapes, nans, axes, dtype, order)
        suite.append(run)

    # partition, argpartition
    funcs = ["partition", "argpartition"]
    for func in funcs:
        if functions is not None and func not in functions:
            continue
        run = {}
        run["name"] = func
        run["statements"] = ["bn_func(a, n, axis)", "sl_func(a, n, axis)"]
        setup = """
            from bottleneck import %s as bn_func
            from bottleneck.slow import %s as sl_func
            if axis is None: n = a.size
            else: n = a.shape[axis] - 1
            n = max(n // 2, 0)
        """ % (
            func,
            func,
        )
        run["setups"] = getsetups(setup, shapes, nans, axes, dtype, order)
        suite.append(run)

    # replace, push
    funcs = ["replace", "push"]
    for func in funcs:
        if functions is not None and func not in functions:
            continue
        run = {}
        run["name"] = func
        if func == "replace":
            run["statements"] = ["bn_func(a, nan, 0)", "slow_func(a, nan, 0)"]
        elif func == "push":
            run["statements"] = ["bn_func(a, 5, axis)", "slow_func(a, 5, axis)"]
        else:
            raise ValueError("Unknow function name")
        setup = """
            from numpy import nan
            from bottleneck import %s as bn_func
            from bottleneck.slow import %s as slow_func
        """ % (
            func,
            func,
        )
        run["setups"] = getsetups(setup, shapes, nans, axes, dtype, order)
        suite.append(run)

    # moving window functions
    funcs = bn.get_functions("move", as_string=True)
    for func in funcs:
        if functions is not None and func not in functions:
            continue
        run = {}
        run["name"] = func
        run["statements"] = ["bn_func(a, w, 1, axis)", "sw_func(a, w, 1, axis)"]
        setup = """
            from bottleneck.slow.move import %s as sw_func
            from bottleneck import %s as bn_func
            w = a.shape[axis] // 5
        """ % (
            func,
            func,
        )
        run["setups"] = getsetups(setup, shapes, nans, axes, dtype, order)
        suite.append(run)

    return suite
