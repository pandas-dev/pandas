from __future__ import annotations

import math
import re
from functools import lru_cache

import numpy as np

from dask.utils import cached_cumsum


@lru_cache(maxsize=512)
def svg(chunks, size=200, **kwargs):
    """Convert chunks from Dask Array into an SVG Image

    Parameters
    ----------
    chunks: tuple
    size: int
        Rough size of the image

    Returns
    -------
    text: An svg string depicting the array as a grid of chunks
    """
    shape = tuple(map(sum, chunks))
    if np.isnan(shape).any():  # don't support unknown sizes
        raise NotImplementedError(
            "Can't generate SVG with unknown chunk sizes.\n\n"
            " A possible solution is with x.compute_chunk_sizes()"
        )
    if not all(shape):
        raise NotImplementedError("Can't generate SVG with 0-length dimensions")
    if len(chunks) == 0:
        raise NotImplementedError("Can't generate SVG with 0 dimensions")
    if len(chunks) == 1:
        return svg_1d(chunks, size=size, **kwargs)
    elif len(chunks) == 2:
        return svg_2d(chunks, size=size, **kwargs)
    elif len(chunks) == 3:
        return svg_3d(chunks, size=size, **kwargs)
    else:
        return svg_nd(chunks, size=size, **kwargs)


text_style = 'font-size="1.0rem" font-weight="100" text-anchor="middle"'


def svg_2d(chunks, offset=(0, 0), skew=(0, 0), size=200, sizes=None):
    shape = tuple(map(sum, chunks))
    sizes = sizes or draw_sizes(shape, size=size)
    y, x = grid_points(chunks, sizes)

    lines, (min_x, max_x, min_y, max_y) = svg_grid(
        x, y, offset=offset, skew=skew, size=size
    )

    header = (
        '<svg width="%d" height="%d" style="stroke:rgb(0,0,0);stroke-width:1" >\n'
        % (max_x + 50, max_y + 50)
    )
    footer = "\n</svg>"

    if shape[0] >= 100:
        rotate = -90
    else:
        rotate = 0

    text = [
        "",
        "  <!-- Text -->",
        '  <text x="%f" y="%f" %s >%d</text>'
        % (max_x / 2, max_y + 20, text_style, shape[1]),
        '  <text x="%f" y="%f" %s transform="rotate(%d,%f,%f)">%d</text>'
        % (max_x + 20, max_y / 2, text_style, rotate, max_x + 20, max_y / 2, shape[0]),
    ]

    return header + "\n".join(lines + text) + footer


def svg_3d(chunks, size=200, sizes=None, offset=(0, 0)):
    shape = tuple(map(sum, chunks))
    sizes = sizes or draw_sizes(shape, size=size)
    x, y, z = grid_points(chunks, sizes)
    ox, oy = offset

    xy, (mnx, mxx, mny, mxy) = svg_grid(
        x / 1.7, y, offset=(ox + 10, oy + 0), skew=(1, 0), size=size
    )

    zx, (_, _, _, max_x) = svg_grid(
        z, x / 1.7, offset=(ox + 10, oy + 0), skew=(0, 1), size=size
    )
    zy, (min_z, max_z, min_y, max_y) = svg_grid(
        z, y, offset=(ox + max_x + 10, oy + max_x), skew=(0, 0), size=size
    )

    header = (
        '<svg width="%d" height="%d" style="stroke:rgb(0,0,0);stroke-width:1" >\n'
        % (max_z + 50, max_y + 50)
    )
    footer = "\n</svg>"

    if shape[1] >= 100:
        rotate = -90
    else:
        rotate = 0

    text = [
        "",
        "  <!-- Text -->",
        '  <text x="%f" y="%f" %s >%d</text>'
        % ((min_z + max_z) / 2, max_y + 20, text_style, shape[2]),
        '  <text x="%f" y="%f" %s transform="rotate(%d,%f,%f)">%d</text>'
        % (
            max_z + 20,
            (min_y + max_y) / 2,
            text_style,
            rotate,
            max_z + 20,
            (min_y + max_y) / 2,
            shape[1],
        ),
        '  <text x="%f" y="%f" %s transform="rotate(45,%f,%f)">%d</text>'
        % (
            (mnx + mxx) / 2 - 10,
            mxy - (mxx - mnx) / 2 + 20,
            text_style,
            (mnx + mxx) / 2 - 10,
            mxy - (mxx - mnx) / 2 + 20,
            shape[0],
        ),
    ]

    return header + "\n".join(xy + zx + zy + text) + footer


def svg_nd(chunks, size=200):
    if len(chunks) % 3 == 1:
        chunks = ((1,),) + chunks
    shape = tuple(map(sum, chunks))
    sizes = draw_sizes(shape, size=size)

    chunks2 = chunks
    sizes2 = sizes
    out = []
    left = 0
    total_height = 0
    while chunks2:
        n = len(chunks2) % 3 or 3
        o = svg(chunks2[:n], sizes=sizes2[:n], offset=(left, 0))
        chunks2 = chunks2[n:]
        sizes2 = sizes2[n:]

        lines = o.split("\n")
        header = lines[0]
        height = float(re.search(r'height="(\d*\.?\d*)"', header).groups()[0])
        total_height = max(total_height, height)
        width = float(re.search(r'width="(\d*\.?\d*)"', header).groups()[0])
        left += width + 10
        o = "\n".join(lines[1:-1])  # remove header and footer

        out.append(o)

    header = (
        '<svg width="%d" height="%d" style="stroke:rgb(0,0,0);stroke-width:1" >\n'
        % (left, total_height)
    )
    footer = "\n</svg>"
    return header + "\n\n".join(out) + footer


def svg_lines(x1, y1, x2, y2, max_n=20):
    """Convert points into lines of text for an SVG plot

    Examples
    --------
    >>> svg_lines([0, 1], [0, 0], [10, 11], [1, 1])  # doctest: +NORMALIZE_WHITESPACE
    ['  <line x1="0" y1="0" x2="10" y2="1" style="stroke-width:2" />',
     '  <line x1="1" y1="0" x2="11" y2="1" style="stroke-width:2" />']
    """
    n = len(x1)

    if n > max_n:
        indices = np.linspace(0, n - 1, max_n, dtype="int")
    else:
        indices = range(n)

    lines = [
        '  <line x1="%d" y1="%d" x2="%d" y2="%d" />' % (x1[i], y1[i], x2[i], y2[i])
        for i in indices
    ]

    lines[0] = lines[0].replace(" /", ' style="stroke-width:2" /')
    lines[-1] = lines[-1].replace(" /", ' style="stroke-width:2" /')
    return lines


def svg_grid(x, y, offset=(0, 0), skew=(0, 0), size=200):
    """Create lines of SVG text that show a grid

    Parameters
    ----------
    x: numpy.ndarray
    y: numpy.ndarray
    offset: tuple
        translational displacement of the grid in SVG coordinates
    skew: tuple
    """
    # Horizontal lines
    x1 = np.zeros_like(y) + offset[0]
    y1 = y + offset[1]
    x2 = np.full_like(y, x[-1]) + offset[0]
    y2 = y + offset[1]

    if skew[0]:
        y2 += x.max() * skew[0]
    if skew[1]:
        x1 += skew[1] * y
        x2 += skew[1] * y

    min_x = min(x1.min(), x2.min())
    min_y = min(y1.min(), y2.min())
    max_x = max(x1.max(), x2.max())
    max_y = max(y1.max(), y2.max())
    max_n = size // 6

    h_lines = ["", "  <!-- Horizontal lines -->"] + svg_lines(x1, y1, x2, y2, max_n)

    # Vertical lines
    x1 = x + offset[0]
    y1 = np.zeros_like(x) + offset[1]
    x2 = x + offset[0]
    y2 = np.full_like(x, y[-1]) + offset[1]

    if skew[0]:
        y1 += skew[0] * x
        y2 += skew[0] * x
    if skew[1]:
        x2 += skew[1] * y.max()

    v_lines = ["", "  <!-- Vertical lines -->"] + svg_lines(x1, y1, x2, y2, max_n)

    color = "ECB172" if len(x) < max_n and len(y) < max_n else "8B4903"
    corners = f"{x1[0]},{y1[0]} {x1[-1]},{y1[-1]} {x2[-1]},{y2[-1]} {x2[0]},{y2[0]}"
    rect = [
        "",
        "  <!-- Colored Rectangle -->",
        f'  <polygon points="{corners}" style="fill:#{color}A0;stroke-width:0"/>',
    ]

    return h_lines + v_lines + rect, (min_x, max_x, min_y, max_y)


def svg_1d(chunks, sizes=None, **kwargs):
    return svg_2d(((1,),) + chunks, **kwargs)


def grid_points(chunks, sizes):
    cumchunks = [np.array(cached_cumsum(c, initial_zero=True)) for c in chunks]
    points = [x * size / x[-1] for x, size in zip(cumchunks, sizes)]
    return points


def draw_sizes(shape, size=200):
    """Get size in pixels for all dimensions"""
    mx = max(shape)
    ratios = [mx / max(0.1, d) for d in shape]
    ratios = [ratio_response(r) for r in ratios]
    return tuple(size / r for r in ratios)


def ratio_response(x):
    """How we display actual size ratios

    Common ratios in sizes span several orders of magnitude,
    which is hard for us to perceive.

    We keep ratios in the 1-3 range accurate, and then apply a logarithm to
    values up until about 100 or so, at which point we stop scaling.
    """
    if x < math.e:
        return x
    elif x <= 100:
        return math.log(x + 12.4)  # f(e) == e
    else:
        return math.log(100 + 12.4)
