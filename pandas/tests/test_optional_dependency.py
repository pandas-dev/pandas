import contextlib
from io import (
    BytesIO,
    StringIO,
)
import sqlite3
import sys
import types

import numpy as np
import pytest

from pandas.compat._optional import (
    VERSIONS,
    import_optional_dependency,
)

import pandas as pd
import pandas._testing as tm


def test_import_optional():
    match = "Import .*notapackage.*pip.*conda.*notapackage"
    with pytest.raises(ImportError, match=match) as exc_info:
        import_optional_dependency("notapackage")
    # The original exception should be there as context:
    assert isinstance(exc_info.value.__context__, ImportError)

    result = import_optional_dependency("notapackage", errors="ignore")
    assert result is None


def test_xlrd_version_fallback():
    pytest.importorskip("xlrd")
    import_optional_dependency("xlrd")


def test_bad_version(monkeypatch):
    name = "fakemodule"
    module = types.ModuleType(name)
    module.__version__ = "0.9.0"
    sys.modules[name] = module
    monkeypatch.setitem(VERSIONS, name, "1.0.0")

    match = "Pandas requires .*1.0.0.* of .fakemodule.*'0.9.0'"
    with pytest.raises(ImportError, match=match):
        import_optional_dependency("fakemodule")

    # Test min_version parameter
    result = import_optional_dependency("fakemodule", min_version="0.8")
    assert result is module

    with tm.assert_produces_warning(UserWarning, match=match):
        result = import_optional_dependency("fakemodule", errors="warn")
    assert result is None

    module.__version__ = "1.0.0"  # exact match is OK
    result = import_optional_dependency("fakemodule")
    assert result is module

    with pytest.raises(ImportError, match="Pandas requires version '1.1.0'"):
        import_optional_dependency("fakemodule", min_version="1.1.0")

    with tm.assert_produces_warning(UserWarning, match="Pandas requires version"):
        result = import_optional_dependency(
            "fakemodule", errors="warn", min_version="1.1.0"
        )
    assert result is None

    result = import_optional_dependency(
        "fakemodule", errors="ignore", min_version="1.1.0"
    )
    assert result is None


def test_submodule(monkeypatch):
    # Create a fake module with a submodule
    name = "fakemodule"
    module = types.ModuleType(name)
    module.__version__ = "0.9.0"
    sys.modules[name] = module
    sub_name = "submodule"
    submodule = types.ModuleType(sub_name)
    setattr(module, sub_name, submodule)
    sys.modules[f"{name}.{sub_name}"] = submodule
    monkeypatch.setitem(VERSIONS, name, "1.0.0")

    match = "Pandas requires .*1.0.0.* of .fakemodule.*'0.9.0'"
    with pytest.raises(ImportError, match=match):
        import_optional_dependency("fakemodule.submodule")

    with tm.assert_produces_warning(UserWarning, match=match):
        result = import_optional_dependency("fakemodule.submodule", errors="warn")
    assert result is None

    module.__version__ = "1.0.0"  # exact match is OK
    result = import_optional_dependency("fakemodule.submodule")
    assert result is submodule


def test_no_version_raises(monkeypatch):
    name = "fakemodule"
    module = types.ModuleType(name)
    sys.modules[name] = module
    monkeypatch.setitem(VERSIONS, name, "1.0.0")

    with pytest.raises(ImportError, match="Can't determine .* fakemodule"):
        import_optional_dependency(name)


# friendlier messages than the bare import error, e.g. GH#19810
CUSTOM_MESSAGES = {
    "matplotlib": "matplotlib is required for plotting",
    "sqlalchemy": "Using a URI string requires 'sqlalchemy'",
}


# These only run where the package is missing, e.g. the CI jobs without optional
# dependencies, which deselect single_cpu tests: don't add module-level skips or marks.
@pytest.mark.parametrize(
    "package, func",
    [
        pytest.param(
            "scipy",
            lambda: pd.Series([1.0, np.nan, 3.0]).interpolate(method="krogh"),
            id="interpolate",
        ),
        pytest.param(
            "scipy",
            lambda: pd.DataFrame({"A": pd.arrays.SparseArray([0, 1])}).sparse.to_coo(),
            id="sparse.to_coo",
        ),
        pytest.param(
            "scipy",
            lambda: pd.Series([1.0, 2.0]).rolling(2, win_type="gaussian").mean(std=1),  # type: ignore[call-arg]
            id="rolling-win_type",
        ),
        pytest.param(
            "scipy",
            lambda: pd.Series([1, 2, 3]).corr(pd.Series([1, 3, 2]), method="kendall"),
            id="Series.corr-kendall",
        ),
        pytest.param(
            "scipy",
            lambda: pd.Series([1, 2, 3]).corr(pd.Series([1, 3, 2]), method="spearman"),
            id="Series.corr-spearman",
        ),
        pytest.param(
            "jinja2", lambda: pd.DataFrame({"a": [1]}).to_latex(), id="to_latex"
        ),
        pytest.param(
            "matplotlib", lambda: pd.DataFrame({"A": [1, 2]}).plot(), id="plot"
        ),
        pytest.param(
            "tabulate", lambda: pd.DataFrame({"a": [1]}).to_markdown(), id="to_markdown"
        ),
        pytest.param(
            "xarray", lambda: pd.DataFrame({"a": [1]}).to_xarray(), id="to_xarray"
        ),
        pytest.param("pyreadstat", lambda: pd.read_spss("data.sav"), id="read_spss"),
        pytest.param(
            "html5lib",
            lambda: pd.read_html(StringIO("<table></table>"), flavor="bs4"),
            id="read_html-bs4",
        ),
        pytest.param(
            "lxml",
            lambda: pd.read_html(StringIO("<table></table>"), flavor="lxml"),
            id="read_html-lxml",
        ),
        pytest.param(
            "openpyxl",
            lambda: pd.read_excel("data.xlsx", engine="openpyxl"),
            id="read_excel-openpyxl",
        ),
        pytest.param(
            "python_calamine",
            lambda: pd.read_excel("data.xlsx", engine="calamine"),
            id="read_excel-calamine",
        ),
        pytest.param(
            "odf",
            lambda: pd.read_excel("data.ods", engine="odf"),
            id="read_excel-odf",
        ),
        pytest.param(
            "xlrd",
            lambda: pd.read_excel("data.xls", engine="xlrd"),
            id="read_excel-xlrd",
        ),
        pytest.param(
            "pyxlsb",
            lambda: pd.read_excel("data.xlsb", engine="pyxlsb"),
            id="read_excel-pyxlsb",
        ),
        pytest.param(
            "openpyxl",
            lambda: pd.DataFrame({"a": [1]}).to_excel(BytesIO(), engine="openpyxl"),
            id="to_excel-openpyxl",
        ),
        pytest.param(
            "xlsxwriter",
            lambda: pd.DataFrame({"a": [1]}).to_excel(BytesIO(), engine="xlsxwriter"),
            id="to_excel-xlsxwriter",
        ),
        pytest.param(
            "odf",
            lambda: pd.DataFrame({"a": [1]}).to_excel(
                BytesIO(),
                engine="odf",  # type: ignore[arg-type]
            ),
            id="to_excel-odf",
        ),
        pytest.param(
            "zstandard",
            lambda: pd.DataFrame({"a": [1]}).to_csv(BytesIO(), compression="zstd"),
            id="to_csv-zstd",
        ),
        pytest.param("pyiceberg", lambda: pd.read_iceberg("t"), id="read_iceberg"),
        pytest.param(
            "pyarrow",
            lambda: pd.read_parquet("data.parquet", engine="pyarrow"),
            id="read_parquet-pyarrow",
        ),
        pytest.param(
            "pyarrow", lambda: pd.read_feather("data.feather"), id="read_feather"
        ),
        pytest.param(
            "pyarrow",
            lambda: pd.DataFrame({"a": [1]}).to_feather(BytesIO()),
            id="to_feather",
        ),
        pytest.param("pyarrow", lambda: pd.read_orc("data.orc"), id="read_orc"),
        pytest.param(
            "numba",
            lambda: (
                pd.Series([1.0, 2.0])
                .rolling(2)
                .apply(lambda x: x.sum(), raw=True, engine="numba")
            ),
            id="rolling.apply-numba",
        ),
        pytest.param(
            "sqlalchemy",
            lambda: pd.read_sql("SELECT 1", "sqlite:///:memory:"),
            id="read_sql-uri",
        ),
    ],
)
def test_missing_optional_dependency_raises(package, func):
    if import_optional_dependency(package, errors="ignore") is not None:
        pytest.skip(f"{package} is installed")
    # messages use the install name, e.g. python-calamine
    match = CUSTOM_MESSAGES.get(package, package.replace("_", "."))
    with pytest.raises(ImportError, match=match):
        func()


def test_style_without_jinja2():
    if import_optional_dependency("jinja2", errors="ignore") is not None:
        pytest.skip("jinja2 is installed")
    # AttributeError so that inspect works without jinja2
    with pytest.raises(AttributeError, match="requires jinja2"):
        pd.DataFrame({"a": [1]}).style


def test_read_sql_unknown_dbapi2_without_sqlalchemy():
    if import_optional_dependency("sqlalchemy", errors="ignore") is not None:
        pytest.skip("sqlalchemy is installed")

    class MockSqliteConnection:
        def __init__(self, *args, **kwargs) -> None:
            self.conn = sqlite3.Connection(*args, **kwargs)

        def __getattr__(self, name):
            return getattr(self.conn, name)

        def close(self):
            self.conn.close()

    with contextlib.closing(MockSqliteConnection(":memory:")) as conn:
        with tm.assert_produces_warning(UserWarning, match="only supports SQLAlchemy"):
            pd.read_sql("SELECT 1", conn)
