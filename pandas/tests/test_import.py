def test_import_error():
    # GH#26367 Verify whether Pandas import succeeds when setting filterwarnings
    import warnings

    warnings.filterwarnings("error")

    try:
        import pandas as pd  # noqa: F401
    except ImportWarning:
        assert False, (
            "Pandas import breaks when settings filterwarnings, see GH#26367. "
            "Please share your issue at "
            "https://github.com/pandas-dev/pandas/issues/26367"
        )
