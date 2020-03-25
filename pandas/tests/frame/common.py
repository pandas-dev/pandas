def _check_mixed_float(df, dtype=None):
    # float16 are most likely to be upcasted to float32
    dtypes = dict(A="float32", B="float32", C="float16", D="float64")
    if isinstance(dtype, str):
        dtypes = {k: dtype for k, v in dtypes.items()}
    elif isinstance(dtype, dict):
        dtypes.update(dtype)
    if dtypes.get("A"):
        assert df.dtypes["A"] == dtypes["A"], (df.dtypes, dtypes)
    if dtypes.get("B"):
        assert df.dtypes["B"] == dtypes["B"], (df.dtypes, dtypes)
    if dtypes.get("C"):
        assert df.dtypes["C"] == dtypes["C"], (df.dtypes, dtypes)
    if dtypes.get("D"):
        assert df.dtypes["D"] == dtypes["D"], (df.dtypes, dtypes)


def _check_mixed_int(df, dtype=None):
    dtypes = dict(A="int32", B="uint64", C="uint8", D="int64")
    if isinstance(dtype, str):
        dtypes = {k: dtype for k, v in dtypes.items()}
    elif isinstance(dtype, dict):
        dtypes.update(dtype)
    if dtypes.get("A"):
        assert df.dtypes["A"] == dtypes["A"], (df.dtypes, dtypes)
    if dtypes.get("B"):
        assert df.dtypes["B"] == dtypes["B"], (df.dtypes, dtypes)
    if dtypes.get("C"):
        assert df.dtypes["C"] == dtypes["C"], (df.dtypes, dtypes)
    if dtypes.get("D"):
        assert df.dtypes["D"] == dtypes["D"], (df.dtypes, dtypes)
