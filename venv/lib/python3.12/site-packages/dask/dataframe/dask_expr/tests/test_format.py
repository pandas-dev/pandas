# flake8: noqa: W291
from __future__ import annotations

from textwrap import dedent

import pytest

from dask.dataframe.dask_expr import from_pandas
from dask.dataframe.dask_expr.tests._util import _backend_library
from dask.dataframe.utils import pyarrow_strings_enabled
from dask.utils import maybe_pluralize

# Set DataFrame backend for this module
pd = _backend_library()


def test_to_string():
    pytest.importorskip("jinja2")
    dtype = "string" if pyarrow_strings_enabled() else "object"
    df = pd.DataFrame(
        {
            "A": [1, 2, 3, 4, 5, 6, 7, 8],
            "B": list("ABCDEFGH"),
            "C": pd.Categorical(list("AAABBBCC")),
        }
    )
    ddf = from_pandas(df, 3)

    exp = dedent(
        f"""\
                           A       B                C
        npartitions=3                                
        0              int64  {dtype}  category[known]
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
      <td>{dtype}</td>
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
</table>"""  # noqa E222, E702
    footer = f"Dask Name: frompandas, {maybe_pluralize(1, 'expression')}"
    exp = f"""<div><strong>Dask DataFrame Structure:</strong></div>
{exp_table}
<div>{footer}</div>"""
    assert ddf.to_html() == exp


def test_series_format():
    s = pd.Series([1, 2, 3, 4, 5, 6, 7, 8], index=list("ABCDEFGH"))
    ds = from_pandas(s, 3)

    exp = dedent(
        """\
    npartitions=3
    A    int64
    D      ...
    G      ...
    H      ..."""
    )
    assert ds.to_string() == exp


def test_series_repr():
    s = pd.Series([1, 2, 3, 4, 5, 6, 7, 8], index=list("ABCDEFGH"))
    dtype = "string" if pyarrow_strings_enabled() else "object"
    ds = from_pandas(s, 3)

    exp = dedent(
        """\
        Dask Series Structure:
        npartitions=3
        A    int64
        D      ...
        G      ...
        H      ...
        Dask Name: frompandas, 1 expression
        Expr=df"""
    )
    assert repr(ds) == exp

    # Not a cheap way to determine if series is empty
    # so does not prefix with "Empty" as we do w/ empty DataFrame
    s = pd.Series([])
    ds = from_pandas(s, 3)

    exp = dedent(
        f"""\
        Dask Series Structure:
        npartitions=3
            {dtype}
               ...
               ...
               ...
        Dask Name: frompandas, 1 expression
        Expr=df"""
    )
    assert repr(ds) == exp


def test_df_repr():
    df = pd.DataFrame({"col1": range(10), "col2": map(float, range(10))})
    ddf = from_pandas(df, 3)

    exp = dedent(
        """\
        Dask DataFrame Structure:
                        col1     col2
        npartitions=3                
        0              int64  float64
        4                ...      ...
        7                ...      ...
        9                ...      ...
        Dask Name: frompandas, 1 expression
        Expr=df"""
    )
    assert repr(ddf) == exp

    df = pd.DataFrame()
    ddf = from_pandas(df, 3)

    exp = dedent(
        """\
        Empty Dask DataFrame Structure:
        npartitions=3
        0              int64  float64
        4                ...      ...
        7                ...      ...
        9                ...      ...
        Dask Name: frompandas, 1 expression
        Expr=df"""
    )


def test_df_to_html():
    pytest.importorskip("jinja2")
    df = pd.DataFrame({"col1": range(10), "col2": map(float, range(10))})
    ddf = from_pandas(df, 3)

    exp = dedent(
        """\
        <div><strong>Dask DataFrame Structure:</strong></div>
        <table border="1" class="dataframe">
          <thead>
            <tr style="text-align: right;">
              <th></th>
              <th>col1</th>
              <th>col2</th>
            </tr>
            <tr>
              <th>npartitions=3</th>
              <th></th>
              <th></th>
            </tr>
          </thead>
          <tbody>
            <tr>
              <th>0</th>
              <td>int64</td>
              <td>float64</td>
            </tr>
            <tr>
              <th>4</th>
              <td>...</td>
              <td>...</td>
            </tr>
            <tr>
              <th>7</th>
              <td>...</td>
              <td>...</td>
            </tr>
            <tr>
              <th>9</th>
              <td>...</td>
              <td>...</td>
            </tr>
          </tbody>
        </table>
        <div>Dask Name: frompandas, 1 expression</div>"""
    )
    assert ddf.to_html() == exp
    assert ddf._repr_html_() == exp  # for jupyter
