# flake8: noqa: W291
from __future__ import annotations

from textwrap import dedent

import numpy as np
import pandas as pd
import pytest

import dask.array as da
import dask.dataframe as dd
from dask.dataframe.utils import get_string_dtype, pyarrow_strings_enabled
from dask.utils import maybe_pluralize

style = """<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
</style>
"""

if dd._dask_expr_enabled():
    pytest.skip("format is completely different", allow_module_level=True)


def _format_string_dtype():
    return "object" if get_string_dtype() is object else "string"


def _format_footer(suffix="", layers=1):
    if pyarrow_strings_enabled():
        return f"Dask Name: to_pyarrow_string{suffix}, {maybe_pluralize(layers + 1, 'graph layer')}"
    return f"Dask Name: from_pandas{suffix}, {maybe_pluralize(layers, 'graph layer')}"


def test_repr():
    df = pd.DataFrame({"x": list(range(100))})
    ddf = dd.from_pandas(df, 3)

    for x in [ddf, ddf.index, ddf.x]:
        assert type(x).__name__ in repr(x)
        assert str(x.npartitions) in repr(x)


def test_repr_meta_mutation():
    # Check that the repr changes when meta changes
    df = pd.DataFrame({"a": range(5), "b": ["a", "b", "c", "d", "e"]})
    ddf = dd.from_pandas(df, npartitions=2)
    s1 = repr(ddf)
    assert repr(ddf) == s1
    ddf.b = ddf.b.astype("category")
    assert repr(ddf) != s1


def test_dataframe_format():
    pytest.importorskip("jinja2")
    df = pd.DataFrame(
        {
            "A": [1, 2, 3, 4, 5, 6, 7, 8],
            "B": list("ABCDEFGH"),
            "C": pd.Categorical(list("AAABBBCC")),
        }
    )
    ddf = dd.from_pandas(df, 3)
    string_dtype = _format_string_dtype()
    footer = _format_footer()
    exp = dedent(
        f"""\
        Dask DataFrame Structure:
                           A       B                C
        npartitions=3                                
        0              int64  {string_dtype}  category[known]
        3                ...     ...              ...
        6                ...     ...              ...
        7                ...     ...              ...
        {footer}"""
    )

    assert repr(ddf) == exp
    assert str(ddf) == exp

    exp = dedent(
        f"""\
                           A       B                C
        npartitions=3                                
        0              int64  {string_dtype}  category[known]
        3                ...     ...              ...
        6                ...     ...              ...
        7                ...     ...              ..."""
    )
    assert ddf.to_string() == exp

    exp_table = f"""<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>A</th>
      <th>B</th>
      <th>C</th>
    </tr>
    <tr>
      <th>npartitions=3</th>
      <th></th>
      <th></th>
      <th></th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>int64</td>
      <td>{string_dtype}</td>
      <td>category[known]</td>
    </tr>
    <tr>
      <th>3</th>
      <td>...</td>
      <td>...</td>
      <td>...</td>
    </tr>
    <tr>
      <th>6</th>
      <td>...</td>
      <td>...</td>
      <td>...</td>
    </tr>
    <tr>
      <th>7</th>
      <td>...</td>
      <td>...</td>
      <td>...</td>
    </tr>
  </tbody>
</table>"""

    exp = f"""<div><strong>Dask DataFrame Structure:</strong></div>
{exp_table}
<div>{footer}</div>"""
    assert ddf.to_html() == exp

    # table is boxed with div and has style
    exp = f"""<div><strong>Dask DataFrame Structure:</strong></div>
<div>
{style}{exp_table}
</div>
<div>{footer}</div>"""
    assert ddf._repr_html_() == exp


def test_dataframe_format_with_index():
    pytest.importorskip("jinja2")
    df = pd.DataFrame(
        {
            "A": [1, 2, 3, 4, 5, 6, 7, 8],
            "B": list("ABCDEFGH"),
            "C": pd.Categorical(list("AAABBBCC")),
        },
        index=list("ABCDEFGH"),
    )
    ddf = dd.from_pandas(df, 3)
    string_dtype = _format_string_dtype()
    footer = _format_footer()
    exp = dedent(
        f"""\
        Dask DataFrame Structure:
                           A       B                C
        npartitions=3                                
        A              int64  {string_dtype}  category[known]
        D                ...     ...              ...
        G                ...     ...              ...
        H                ...     ...              ...
        {footer}"""
    )
    assert repr(ddf) == exp
    assert str(ddf) == exp

    exp_table = f"""<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>A</th>
      <th>B</th>
      <th>C</th>
    </tr>
    <tr>
      <th>npartitions=3</th>
      <th></th>
      <th></th>
      <th></th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>A</th>
      <td>int64</td>
      <td>{string_dtype}</td>
      <td>category[known]</td>
    </tr>
    <tr>
      <th>D</th>
      <td>...</td>
      <td>...</td>
      <td>...</td>
    </tr>
    <tr>
      <th>G</th>
      <td>...</td>
      <td>...</td>
      <td>...</td>
    </tr>
    <tr>
      <th>H</th>
      <td>...</td>
      <td>...</td>
      <td>...</td>
    </tr>
  </tbody>
</table>"""

    exp = f"""<div><strong>Dask DataFrame Structure:</strong></div>
{exp_table}
<div>{footer}</div>"""
    assert ddf.to_html() == exp

    # table is boxed with div and has style
    exp = f"""<div><strong>Dask DataFrame Structure:</strong></div>
<div>
{style}{exp_table}
</div>
<div>{footer}</div>"""
    assert ddf._repr_html_() == exp


def test_dataframe_format_unknown_divisions():
    pytest.importorskip("jinja2")
    df = pd.DataFrame(
        {
            "A": [1, 2, 3, 4, 5, 6, 7, 8],
            "B": list("ABCDEFGH"),
            "C": pd.Categorical(list("AAABBBCC")),
        }
    )
    ddf = dd.from_pandas(df, 3)
    ddf = ddf.clear_divisions()
    assert not ddf.known_divisions

    string_dtype = _format_string_dtype()
    footer = _format_footer()

    exp = dedent(
        f"""\
        Dask DataFrame Structure:
                           A       B                C
        npartitions=3                                
                       int64  {string_dtype}  category[known]
                         ...     ...              ...
                         ...     ...              ...
                         ...     ...              ...
        {footer}"""
    )
    assert repr(ddf) == exp
    assert str(ddf) == exp

    exp = dedent(
        f"""\
                           A       B                C
        npartitions=3                                
                       int64  {string_dtype}  category[known]
                         ...     ...              ...
                         ...     ...              ...
                         ...     ...              ..."""
    )
    assert ddf.to_string() == exp

    exp_table = f"""<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>A</th>
      <th>B</th>
      <th>C</th>
    </tr>
    <tr>
      <th>npartitions=3</th>
      <th></th>
      <th></th>
      <th></th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th></th>
      <td>int64</td>
      <td>{string_dtype}</td>
      <td>category[known]</td>
    </tr>
    <tr>
      <th></th>
      <td>...</td>
      <td>...</td>
      <td>...</td>
    </tr>
    <tr>
      <th></th>
      <td>...</td>
      <td>...</td>
      <td>...</td>
    </tr>
    <tr>
      <th></th>
      <td>...</td>
      <td>...</td>
      <td>...</td>
    </tr>
  </tbody>
</table>"""

    exp = f"""<div><strong>Dask DataFrame Structure:</strong></div>
{exp_table}
<div>{footer}</div>"""
    assert ddf.to_html() == exp

    # table is boxed with div and has style
    exp = f"""<div><strong>Dask DataFrame Structure:</strong></div>
<div>
{style}{exp_table}
</div>
<div>{footer}</div>"""
    assert ddf._repr_html_() == exp


def test_dataframe_format_long():
    pytest.importorskip("jinja2")
    df = pd.DataFrame(
        {
            "A": [1, 2, 3, 4, 5, 6, 7, 8] * 10,
            "B": list("ABCDEFGH") * 10,
            "C": pd.Categorical(list("AAABBBCC") * 10),
        }
    )
    string_dtype = _format_string_dtype()
    footer = _format_footer()
    ddf = dd.from_pandas(df, 10)
    exp = dedent(
        f"""\
        Dask DataFrame Structure:
                            A       B                C
        npartitions=10                                
        0               int64  {string_dtype}  category[known]
        8                 ...     ...              ...
        ...               ...     ...              ...
        72                ...     ...              ...
        79                ...     ...              ...
        {footer}"""
    )
    assert repr(ddf) == exp
    assert str(ddf) == exp

    exp = dedent(
        f"""\
                            A       B                C
        npartitions=10                                
        0               int64  {string_dtype}  category[known]
        8                 ...     ...              ...
        ...               ...     ...              ...
        72                ...     ...              ...
        79                ...     ...              ..."""
    )
    assert ddf.to_string() == exp

    exp_table = f"""<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>A</th>
      <th>B</th>
      <th>C</th>
    </tr>
    <tr>
      <th>npartitions=10</th>
      <th></th>
      <th></th>
      <th></th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>int64</td>
      <td>{string_dtype}</td>
      <td>category[known]</td>
    </tr>
    <tr>
      <th>8</th>
      <td>...</td>
      <td>...</td>
      <td>...</td>
    </tr>
    <tr>
      <th>...</th>
      <td>...</td>
      <td>...</td>
      <td>...</td>
    </tr>
    <tr>
      <th>72</th>
      <td>...</td>
      <td>...</td>
      <td>...</td>
    </tr>
    <tr>
      <th>79</th>
      <td>...</td>
      <td>...</td>
      <td>...</td>
    </tr>
  </tbody>
</table>"""

    exp = f"""<div><strong>Dask DataFrame Structure:</strong></div>
{exp_table}
<div>{footer}</div>"""
    assert ddf.to_html() == exp

    # table is boxed with div
    exp = f"""<div><strong>Dask DataFrame Structure:</strong></div>
<div>
{style}{exp_table}
</div>
<div>{footer}</div>"""
    assert ddf._repr_html_() == exp


def test_series_format():
    s = pd.Series([1, 2, 3, 4, 5, 6, 7, 8], index=list("ABCDEFGH"))
    ds = dd.from_pandas(s, 3)
    footer = _format_footer()
    exp = dedent(
        f"""\
    Dask Series Structure:
    npartitions=3
    A    int64
    D      ...
    G      ...
    H      ...
    dtype: int64
    {footer}"""
    )
    assert repr(ds) == exp
    assert str(ds) == exp

    exp = dedent(
        """\
    npartitions=3
    A    int64
    D      ...
    G      ...
    H      ..."""
    )
    assert ds.to_string() == exp

    s = pd.Series([1, 2, 3, 4, 5, 6, 7, 8], index=list("ABCDEFGH"), name="XXX")
    ds = dd.from_pandas(s, 3)
    exp = dedent(
        f"""\
    Dask Series Structure:
    npartitions=3
    A    int64
    D      ...
    G      ...
    H      ...
    Name: XXX, dtype: int64
    {footer}"""
    )
    assert repr(ds) == exp
    assert str(ds) == exp


def test_series_format_long():
    s = pd.Series([1, 2, 3, 4, 5, 6, 7, 8, 9, 10] * 10, index=list("ABCDEFGHIJ") * 10)
    ds = dd.from_pandas(s, 10)
    footer = _format_footer()
    exp = dedent(
        f"""\
        Dask Series Structure:
        npartitions=10
        A    int64
        B      ...
             ...  
        J      ...
        J      ...
        dtype: int64
        {footer}"""
    )
    assert repr(ds) == exp
    assert str(ds) == exp
    exp = dedent(
        """\
    npartitions=10
    A    int64
    B      ...
         ...  
    J      ...
    J      ..."""
    )
    assert ds.to_string() == exp


def test_index_format():
    s = pd.Series([1, 2, 3, 4, 5, 6, 7, 8], index=list("ABCDEFGH"))
    ds = dd.from_pandas(s, 3)
    string_dtype = _format_string_dtype()
    footer = _format_footer("-index", 2)
    exp = dedent(
        f"""\
    Dask Index Structure:
    npartitions=3
    A    {string_dtype}
    D       ...
    G       ...
    H       ...
    dtype: {string_dtype}
    {footer}"""
    )
    assert repr(ds.index) == exp
    assert str(ds.index) == exp

    s = pd.Series(
        [1, 2, 3, 4, 5, 6, 7, 8],
        index=pd.CategoricalIndex([1, 2, 3, 4, 5, 6, 7, 8], name="YYY"),
    )
    ds = dd.from_pandas(s, 3)
    exp = dedent(
        """\
    Dask Index Structure:
    npartitions=3
    1    category[known]
    4                ...
    7                ...
    8                ...
    Name: YYY, dtype: category
    Dask Name: from_pandas-index, 2 graph layers"""
    )
    assert repr(ds.index) == exp
    assert str(ds.index) == exp


def test_categorical_format():
    s = pd.Series(["a", "b", "c"]).astype("category")
    known = dd.from_pandas(s, npartitions=1)
    unknown = known.cat.as_unknown()
    exp = (
        "Dask Series Structure:\n"
        "npartitions=1\n"
        "0    category[known]\n"
        "2                ...\n"
        "dtype: category\n"
        "Dask Name: from_pandas, 1 graph layer"
    )
    assert repr(known) == exp
    exp = (
        "Dask Series Structure:\n"
        "npartitions=1\n"
        "0    category[unknown]\n"
        "2                  ...\n"
        "dtype: category\n"
        "Dask Name: from_pandas, 1 graph layer"
    )
    assert repr(unknown) == exp


def test_duplicate_columns_repr():
    arr = da.from_array(np.arange(10).reshape(5, 2), chunks=(5, 2))
    frame = dd.from_dask_array(arr, columns=["a", "a"])
    repr(frame)


def test_empty_repr():
    pytest.importorskip("jinja2")
    df = pd.DataFrame()
    ddf = dd.from_pandas(df, npartitions=1)
    exp = (
        "Empty Dask DataFrame Structure:\n"
        "Columns: []\n"
        "Divisions: [, ]\n"
        "Dask Name: from_pandas, 1 graph layer"
    )
    assert repr(ddf) == exp
    exp_table = """<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
    </tr>
    <tr>
      <th>npartitions=1</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th></th>
    </tr>
    <tr>
      <th></th>
    </tr>
  </tbody>
</table>"""
    exp = """<div><strong>Dask DataFrame Structure:</strong></div>
<div>
{style}{exp_table}
</div>
<div>Dask Name: from_pandas, 1 graph layer</div>""".format(
        style=style, exp_table=exp_table
    )
    assert ddf._repr_html_() == exp
