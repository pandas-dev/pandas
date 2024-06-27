from __future__ import annotations

import xml.etree.ElementTree

import pytest

import dask.array as da
from dask.array.svg import draw_sizes


def parses(text):
    cleaned = text.replace("&rarr;", "")  # xml doesn't like righarrow character
    assert xml.etree.ElementTree.fromstring(cleaned) is not None  # parses cleanly


def test_basic():
    parses(da.ones(10).to_svg())
    parses(da.ones((10, 10)).to_svg())
    parses(da.ones((10, 10, 10)).to_svg())
    parses(da.ones((10, 10, 10, 10)).to_svg())
    parses(da.ones((10, 10, 10, 10, 10)).to_svg())
    parses(da.ones((10, 10, 10, 10, 10, 10)).to_svg())
    parses(da.ones((10, 10, 10, 10, 10, 10, 10)).to_svg())


def test_repr_html():
    pytest.importorskip("jinja2")
    assert da.ones([])._repr_html_()
    assert da.ones(10)[:0]._repr_html_()
    assert da.ones(10)._repr_html_()
    assert da.ones((10, 10))._repr_html_()
    assert da.ones((10, 10, 10))._repr_html_()
    assert da.ones((10, 10, 10, 10))._repr_html_()


def test_errors():
    # empty arrays
    with pytest.raises(NotImplementedError) as excpt:
        da.ones([]).to_svg()
    assert "0 dimensions" in str(excpt.value)

    # Scalars
    with pytest.raises(NotImplementedError) as excpt:
        da.asarray(1).to_svg()
    assert "0 dimensions" in str(excpt.value)

    # 0-length dims arrays
    with pytest.raises(NotImplementedError) as excpt:
        da.ones(10)[:0].to_svg()
    assert "0-length dimensions" in str(excpt.value)

    # unknown chunk sizes
    with pytest.raises(NotImplementedError) as excpt:
        x = da.ones(10)
        x = x[x > 5]
        x.to_svg()
    assert "unknown chunk sizes" in str(excpt.value)


def test_repr_html_size_units():
    pytest.importorskip("jinja2")
    x = da.ones((10000, 5000))
    x = da.ones((3000, 10000), chunks=(1000, 1000))
    text = x._repr_html_()

    assert "MB" in text or "MiB" in text
    assert str(x.shape) in text
    assert str(x.dtype) in text

    parses(text)

    x = da.ones((3000, 10000, 50), chunks=(1000, 1000, 10))
    parses(x._repr_html_())


def test_draw_sizes():
    assert draw_sizes((10, 10), size=100) == (100, 100)  # respect symmetry
    assert draw_sizes((10, 10), size=200) == (200, 200)  # respect size keyword
    assert draw_sizes((10, 5), size=100) == (100, 50)  # respect small ratios

    a, b, c = draw_sizes((1000, 100, 10))
    assert a > b
    assert b > c
    assert a < b * 5
    assert b < c * 5


def test_too_many_lines_fills_sides_darker():
    data = da.ones((16000, 2400, 3600), chunks=(1, 2400, 3600))
    text = data.to_svg()
    assert "8B4903" in text
    assert text.count("\n") < 300


def test_3d():
    text = da.ones((10, 10, 10, 10, 10)).to_svg()
    assert text.count("<svg") == 1
