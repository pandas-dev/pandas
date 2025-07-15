from __future__ import annotations

import importlib
import warnings

# The "array.query-planning" config can only be processed once
ARRAY_EXPR_ENABLED: bool | None = None


def _array_expr_enabled() -> bool:
    import dask

    global ARRAY_EXPR_ENABLED

    use_array_expr = dask.config.get("array.query-planning")

    if ARRAY_EXPR_ENABLED is not None:
        if (use_array_expr is True and ARRAY_EXPR_ENABLED is False) or (
            use_array_expr is False and ARRAY_EXPR_ENABLED is True
        ):
            warnings.warn(
                "The 'array.query-planning' config is now set to "
                f"{use_array_expr}, but query planning is already "
                f"{'enabled' if ARRAY_EXPR_ENABLED else 'disabled'}. "
                "The query-planning config can only be changed before "
                "`dask.array` is first imported!"
            )
        return ARRAY_EXPR_ENABLED  # type: ignore[return-value]

    return bool(use_array_expr if use_array_expr is not None else False)


def array_expr_enabled():
    # Need a public variant for downstream libraries to check
    return _array_expr_enabled()


__all__ = [
    "bool",
    "complex64",
    "complex128",
    "e",
    "euler_gamma",
    "float32",
    "float64",
    "inf",
    "int8",
    "int16",
    "int32",
    "int64",
    "nan",
    "newaxis",
    "pi",
    "uint8",
    "uint16",
    "uint32",
    "uint64",
    "backends",
    "fft",
    "lib",
    "linalg",
    "ma",
    "overlap",
    "random",
    "array_expr_enabled",
    "shuffle",
    "atop",
    "blockwise",
    "register_chunk_type",
    "Array",
    "PerformanceWarning",
    "asanyarray",
    "asarray",
    "block",
    "broadcast_arrays",
    "broadcast_to",
    "concatenate",
    "from_array",
    "from_delayed",
    "from_npy_stack",
    "from_zarr",
    "map_blocks",
    "stack",
    "store",
    "to_hdf5",
    "to_npy_stack",
    "to_zarr",
    "unify_chunks",
    "arange",
    "diag",
    "diagonal",
    "empty_like",
    "eye",
    "fromfunction",
    "full_like",
    "indices",
    "linspace",
    "meshgrid",
    "ones_like",
    "pad",
    "repeat",
    "tile",
    "tri",
    "zeros_like",
    "apply_gufunc",
    "as_gufunc",
    "gufunc",
    "moveaxis",
    "rollaxis",
    "optimize",
    "map_overlap",
    "push",
    "nanpercentile",
    "percentile",
    "rechunk",
    "all",
    "any",
    "argmax",
    "argmin",
    "argtopk",
    "cumprod",
    "cumsum",
    "max",
    "mean",
    "median",
    "min",
    "moment",
    "nanargmax",
    "nanargmin",
    "nancumprod",
    "nancumsum",
    "nanmax",
    "nanmean",
    "nanmedian",
    "nanmin",
    "nanprod",
    "nanquantile",
    "nanstd",
    "nansum",
    "nanvar",
    "prod",
    "quantile",
    "reduction",
    "std",
    "sum",
    "topk",
    "trace",
    "var",
    "reshape",
    "reshape_blockwise",
    "allclose",
    "append",
    "apply_along_axis",
    "apply_over_axes",
    "argwhere",
    "around",
    "array",
    "atleast_1d",
    "atleast_2d",
    "atleast_3d",
    "average",
    "bincount",
    "choose",
    "coarsen",
    "compress",
    "corrcoef",
    "count_nonzero",
    "cov",
    "delete",
    "diff",
    "digitize",
    "dot",
    "dstack",
    "ediff1d",
    "einsum",
    "expand_dims",
    "extract",
    "flatnonzero",
    "flip",
    "fliplr",
    "flipud",
    "gradient",
    "histogram",
    "histogram2d",
    "histogramdd",
    "hstack",
    "insert",
    "isclose",
    "isin",
    "isnull",
    "matmul",
    "ndim",
    "nonzero",
    "notnull",
    "outer",
    "piecewise",
    "ptp",
    "ravel",
    "ravel_multi_index",
    "result_type",
    "roll",
    "rot90",
    "round",
    "searchsorted",
    "select",
    "shape",
    "squeeze",
    "swapaxes",
    "take",
    "tensordot",
    "transpose",
    "tril",
    "tril_indices",
    "tril_indices_from",
    "triu",
    "triu_indices",
    "triu_indices_from",
    "union1d",
    "unique",
    "unravel_index",
    "vdot",
    "vstack",
    "where",
    "from_tiledb",
    "to_tiledb",
    "abs",
    "absolute",
    "add",
    "angle",
    "arccos",
    "arccosh",
    "arcsin",
    "arcsinh",
    "arctan",
    "arctan2",
    "arctanh",
    "bitwise_and",
    "bitwise_not",
    "bitwise_or",
    "bitwise_xor",
    "cbrt",
    "ceil",
    "clip",
    "conj",
    "copysign",
    "cos",
    "cosh",
    "deg2rad",
    "degrees",
    "divide",
    "divmod",
    "equal",
    "exp",
    "exp2",
    "expm1",
    "fabs",
    "fix",
    "float_power",
    "floor",
    "floor_divide",
    "fmax",
    "fmin",
    "fmod",
    "frexp",
    "frompyfunc",
    "greater",
    "greater_equal",
    "hypot",
    "i0",
    "imag",
    "invert",
    "iscomplex",
    "isfinite",
    "isinf",
    "isnan",
    "isneginf",
    "isposinf",
    "isreal",
    "ldexp",
    "left_shift",
    "less",
    "less_equal",
    "log",
    "log1p",
    "log2",
    "log10",
    "logaddexp",
    "logaddexp2",
    "logical_and",
    "logical_not",
    "logical_or",
    "logical_xor",
    "maximum",
    "minimum",
    "mod",
    "modf",
    "multiply",
    "nan_to_num",
    "negative",
    "nextafter",
    "not_equal",
    "positive",
    "power",
    "rad2deg",
    "radians",
    "real",
    "reciprocal",
    "remainder",
    "right_shift",
    "rint",
    "sign",
    "signbit",
    "sin",
    "sinc",
    "sinh",
    "spacing",
    "sqrt",
    "square",
    "subtract",
    "tan",
    "tanh",
    "true_divide",
    "trunc",
    "assert_eq",
    "empty",
    "full",
    "ones",
    "zeros",
    "compute",
]

try:
    from numpy import bool_ as bool
    from numpy import (
        complex64,
        complex128,
        e,
        euler_gamma,
        float32,
        float64,
        inf,
        int8,
        int16,
        int32,
        int64,
        nan,
        newaxis,
        pi,
        uint8,
        uint16,
        uint32,
        uint64,
    )

    from dask.array import backends, fft, lib, linalg, ma, overlap, random
    from dask.array._shuffle import shuffle
    from dask.array.blockwise import atop, blockwise
    from dask.array.chunk_types import register_chunk_type
    from dask.array.core import (
        Array,
        PerformanceWarning,
        asanyarray,
        asarray,
        block,
        broadcast_arrays,
        broadcast_to,
        concatenate,
        from_array,
        from_delayed,
        from_npy_stack,
        from_zarr,
        map_blocks,
        stack,
        store,
        to_hdf5,
        to_npy_stack,
        to_zarr,
        unify_chunks,
    )
    from dask.array.creation import (
        arange,
        diag,
        diagonal,
        empty_like,
        eye,
        fromfunction,
        full_like,
        indices,
        linspace,
        meshgrid,
        ones_like,
        pad,
        repeat,
        tile,
        tri,
        zeros_like,
    )
    from dask.array.gufunc import apply_gufunc, as_gufunc, gufunc
    from dask.array.numpy_compat import moveaxis, rollaxis
    from dask.array.optimization import optimize
    from dask.array.overlap import map_overlap, push
    from dask.array.percentile import nanpercentile, percentile
    from dask.array.rechunk import rechunk
    from dask.array.reductions import (
        all,
        any,
        argmax,
        argmin,
        argtopk,
        cumprod,
        cumsum,
        max,
        mean,
        median,
        min,
        moment,
        nanargmax,
        nanargmin,
        nancumprod,
        nancumsum,
        nanmax,
        nanmean,
        nanmedian,
        nanmin,
        nanprod,
        nanquantile,
        nanstd,
        nansum,
        nanvar,
        prod,
        quantile,
        reduction,
        std,
        sum,
        topk,
        trace,
        var,
    )
    from dask.array.reshape import reshape, reshape_blockwise
    from dask.array.routines import (
        allclose,
        append,
        apply_along_axis,
        apply_over_axes,
        argwhere,
        around,
        array,
        atleast_1d,
        atleast_2d,
        atleast_3d,
        average,
        bincount,
        choose,
        coarsen,
        compress,
        corrcoef,
        count_nonzero,
        cov,
        delete,
        diff,
        digitize,
        dot,
        dstack,
        ediff1d,
        einsum,
        expand_dims,
        extract,
        flatnonzero,
        flip,
        fliplr,
        flipud,
        gradient,
        histogram,
        histogram2d,
        histogramdd,
        hstack,
        insert,
        isclose,
        isin,
        isnull,
        matmul,
        ndim,
        nonzero,
        notnull,
        outer,
        piecewise,
        ptp,
        ravel,
        ravel_multi_index,
        result_type,
        roll,
        rot90,
        round,
        searchsorted,
        select,
        shape,
        squeeze,
        swapaxes,
        take,
        tensordot,
        transpose,
        tril,
        tril_indices,
        tril_indices_from,
        triu,
        triu_indices,
        triu_indices_from,
        union1d,
        unique,
        unravel_index,
        vdot,
        vstack,
        where,
    )
    from dask.array.tiledb_io import from_tiledb, to_tiledb
    from dask.array.ufunc import (
        abs,
        absolute,
        add,
        angle,
        arccos,
        arccosh,
        arcsin,
        arcsinh,
        arctan,
        arctan2,
        arctanh,
        bitwise_and,
        bitwise_not,
        bitwise_or,
        bitwise_xor,
        cbrt,
        ceil,
        clip,
        conj,
        copysign,
        cos,
        cosh,
        deg2rad,
        degrees,
        divide,
        divmod,
        equal,
        exp,
        exp2,
        expm1,
        fabs,
        fix,
        float_power,
        floor,
        floor_divide,
        fmax,
        fmin,
        fmod,
        frexp,
        frompyfunc,
        greater,
        greater_equal,
        hypot,
        i0,
        imag,
        invert,
        iscomplex,
        isfinite,
        isinf,
        isnan,
        isneginf,
        isposinf,
        isreal,
        ldexp,
        left_shift,
        less,
        less_equal,
        log,
        log1p,
        log2,
        log10,
        logaddexp,
        logaddexp2,
        logical_and,
        logical_not,
        logical_or,
        logical_xor,
        maximum,
        minimum,
        mod,
        modf,
        multiply,
        nan_to_num,
        negative,
        nextafter,
        not_equal,
        positive,
        power,
        rad2deg,
        radians,
        real,
        reciprocal,
        remainder,
        right_shift,
        rint,
        sign,
        signbit,
        sin,
        sinc,
        sinh,
        spacing,
        sqrt,
        square,
        subtract,
        tan,
        tanh,
        true_divide,
        trunc,
    )
    from dask.array.utils import assert_eq
    from dask.array.wrap import empty, full, ones, zeros
    from dask.base import compute

    if _array_expr_enabled():
        import dask.array._array_expr as da

        da = importlib.reload(da)

except ImportError as e:
    msg = (
        "Dask array requirements are not installed.\n\n"
        "Please either conda or pip install as follows:\n\n"
        "  conda install dask                 # either conda install\n"
        '  python -m pip install "dask[array]" --upgrade  # or python -m pip install'
    )
    raise ImportError(str(e) + "\n\n" + msg) from e


if _array_expr_enabled():

    def raise_not_implemented_error(attr_name):
        def inner_func(*args, **kwargs):
            raise NotImplementedError(
                f"Function {attr_name} is not implemented for dask-expr."
            )

        return inner_func

    try:
        from dask.array._array_expr import Array  # type: ignore
        from dask.array._array_expr import _overlap as overlap  # type: ignore
        from dask.array._array_expr import (  # type: ignore
            abs,
            absolute,
            add,
            angle,
            apply_gufunc,
            arange,
            arccos,
            arccosh,
            arcsin,
            arcsinh,
            arctan,
            arctan2,
            arctanh,
            array,
            as_gufunc,
            asanyarray,
            asarray,
            bitwise_and,
            bitwise_not,
            bitwise_or,
            bitwise_xor,
            blockwise,
            cbrt,
            ceil,
            clip,
            concatenate,
            conj,
            copysign,
            cos,
            cosh,
            deg2rad,
            degrees,
            divide,
            divmod,
            elemwise,
            empty,
            empty_like,
            equal,
            exp,
            exp2,
            expm1,
            fabs,
            fix,
            float_power,
            floor,
            floor_divide,
            fmax,
            fmin,
            fmod,
            frexp,
            from_array,
            frompyfunc,
            full,
            full_like,
            greater,
            greater_equal,
            gufunc,
            hypot,
            i0,
            imag,
            invert,
            iscomplex,
            isfinite,
            isinf,
            isnan,
            isneginf,
            isposinf,
            isreal,
            ldexp,
            left_shift,
            less,
            less_equal,
            linspace,
            log,
            log1p,
            log2,
            log10,
            logaddexp,
            logaddexp2,
            logical_and,
            logical_not,
            logical_or,
            logical_xor,
            map_blocks,
            map_overlap,
            maximum,
            minimum,
            mod,
            modf,
            multiply,
            nan_to_num,
            negative,
            nextafter,
            not_equal,
            ones,
            ones_like,
            positive,
            power,
            rad2deg,
            radians,
            random,
            real,
            rechunk,
            reciprocal,
            reduction,
            remainder,
            repeat,
            right_shift,
            rint,
            sign,
            signbit,
            sin,
            sinc,
            sinh,
            spacing,
            sqrt,
            square,
            stack,
            subtract,
            tan,
            tanh,
            true_divide,
            trunc,
            zeros,
            zeros_like,
        )
        from dask.array.reductions import (
            all,
            any,
            max,
            mean,
            min,
            moment,
            nanmax,
            nanmean,
            nanmin,
            nanprod,
            nanstd,
            nansum,
            nanvar,
            prod,
            reduction,
            std,
            sum,
            var,
        )

        backends = raise_not_implemented_error("backends")
        fft = raise_not_implemented_error("fft")
        lib = raise_not_implemented_error("lib")
        linalg = raise_not_implemented_error("linalg")
        ma = raise_not_implemented_error("ma")
        atop = raise_not_implemented_error("atop")
        register_chunk_type = raise_not_implemented_error("register_chunk_type")
        block = raise_not_implemented_error("block")
        broadcast_arrays = raise_not_implemented_error("broadcast_arrays")
        broadcast_to = raise_not_implemented_error("broadcast_to")
        from_delayed = raise_not_implemented_error("from_delayed")
        from_npy_stack = raise_not_implemented_error("from_npy_stack")
        from_zarr = raise_not_implemented_error("from_zarr")
        store = raise_not_implemented_error("store")
        to_hdf5 = raise_not_implemented_error("to_hdf5")
        to_npy_stack = raise_not_implemented_error("to_npy_stack")
        to_zarr = raise_not_implemented_error("to_zarr")
        unify_chunks = raise_not_implemented_error("unify_chunks")
        diag = raise_not_implemented_error("diag")
        diagonal = raise_not_implemented_error("diagonal")
        eye = raise_not_implemented_error("eye")
        fromfunction = raise_not_implemented_error("fromfunction")
        indices = raise_not_implemented_error("indices")
        meshgrid = raise_not_implemented_error("meshgrid")
        pad = raise_not_implemented_error("pad")
        tile = raise_not_implemented_error("tile")
        tri = raise_not_implemented_error("tri")
        moveaxis = raise_not_implemented_error("moveaxis")
        rollaxis = raise_not_implemented_error("rollaxis")
        optimize = raise_not_implemented_error("optimize")
        percentile = raise_not_implemented_error("percentile")
        argmax = raise_not_implemented_error("argmax")
        argmin = raise_not_implemented_error("argmin")
        argtopk = raise_not_implemented_error("argtopk")
        cumprod = raise_not_implemented_error("cumprod")
        cumsum = raise_not_implemented_error("cumsum")
        median = raise_not_implemented_error("median")
        nanargmax = raise_not_implemented_error("nanargmax")
        nanargmin = raise_not_implemented_error("nanargmin")
        nancumprod = raise_not_implemented_error("nancumprod")
        nancumsum = raise_not_implemented_error("nancumsum")
        nanmedian = raise_not_implemented_error("nanmedian")
        topk = raise_not_implemented_error("topk")
        trace = raise_not_implemented_error("trace")
        reshape = raise_not_implemented_error("reshape")
        allclose = raise_not_implemented_error("allclose")
        append = raise_not_implemented_error("append")
        apply_along_axis = raise_not_implemented_error("apply_along_axis")
        apply_over_axes = raise_not_implemented_error("apply_over_axes")
        argwhere = raise_not_implemented_error("argwhere")
        around = raise_not_implemented_error("around")
        atleast_1d = raise_not_implemented_error("atleast_1d")
        atleast_2d = raise_not_implemented_error("atleast_2d")
        atleast_3d = raise_not_implemented_error("atleast_3d")
        average = raise_not_implemented_error("average")
        bincount = raise_not_implemented_error("bincount")
        choose = raise_not_implemented_error("choose")
        coarsen = raise_not_implemented_error("coarsen")
        compress = raise_not_implemented_error("compress")
        corrcoef = raise_not_implemented_error("corrcoef")
        count_nonzero = raise_not_implemented_error("count_nonzero")
        cov = raise_not_implemented_error("cov")
        delete = raise_not_implemented_error("delete")
        diff = raise_not_implemented_error("diff")
        digitize = raise_not_implemented_error("digitize")
        dot = raise_not_implemented_error("dot")
        dstack = raise_not_implemented_error("dstack")
        ediff1d = raise_not_implemented_error("ediff1d")
        einsum = raise_not_implemented_error("einsum")
        expand_dims = raise_not_implemented_error("expand_dims")
        extract = raise_not_implemented_error("extract")
        flatnonzero = raise_not_implemented_error("flatnonzero")
        flip = raise_not_implemented_error("flip")
        fliplr = raise_not_implemented_error("fliplr")
        flipud = raise_not_implemented_error("flipud")
        gradient = raise_not_implemented_error("gradient")
        histogram = raise_not_implemented_error("histogram")
        histogram2d = raise_not_implemented_error("histogram2d")
        histogramdd = raise_not_implemented_error("histogramdd")
        hstack = raise_not_implemented_error("hstack")
        insert = raise_not_implemented_error("insert")
        isclose = raise_not_implemented_error("isclose")
        isin = raise_not_implemented_error("isin")
        isnull = raise_not_implemented_error("isnull")
        matmul = raise_not_implemented_error("matmul")
        ndim = raise_not_implemented_error("ndim")
        nonzero = raise_not_implemented_error("nonzero")
        notnull = raise_not_implemented_error("notnull")
        outer = raise_not_implemented_error("outer")
        piecewise = raise_not_implemented_error("piecewise")
        ptp = raise_not_implemented_error("ptp")
        ravel = raise_not_implemented_error("ravel")
        ravel_multi_index = raise_not_implemented_error("ravel_multi_index")
        result_type = raise_not_implemented_error("result_type")
        roll = raise_not_implemented_error("roll")
        rot90 = raise_not_implemented_error("rot90")
        round = raise_not_implemented_error("round")
        searchsorted = raise_not_implemented_error("searchsorted")
        select = raise_not_implemented_error("select")
        shape = raise_not_implemented_error("shape")
        squeeze = raise_not_implemented_error("squeeze")
        swapaxes = raise_not_implemented_error("swapaxes")
        take = raise_not_implemented_error("take")
        tensordot = raise_not_implemented_error("tensordot")
        transpose = raise_not_implemented_error("transpose")
        tril = raise_not_implemented_error("tril")
        tril_indices = raise_not_implemented_error("tril_indices")
        tril_indices_from = raise_not_implemented_error("tril_indices_from")
        triu = raise_not_implemented_error("triu")
        triu_indices = raise_not_implemented_error("triu_indices")
        triu_indices_from = raise_not_implemented_error("triu_indices_from")
        union1d = raise_not_implemented_error("union1d")
        unique = raise_not_implemented_error("unique")
        unravel_index = raise_not_implemented_error("unravel_index")
        vdot = raise_not_implemented_error("vdot")
        vstack = raise_not_implemented_error("vstack")
        where = raise_not_implemented_error("where")
        from_tiledb = raise_not_implemented_error("from_tiledb")
        to_tiledb = raise_not_implemented_error("to_tiledb")

        from dask.array.utils import assert_eq
        from dask.base import compute

    except ImportError:
        import dask.array as da  # type: ignore

        da = importlib.reload(da)
