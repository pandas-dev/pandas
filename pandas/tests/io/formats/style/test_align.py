import numpy as np
import pytest

from pandas import DataFrame

pytest.importorskip("jinja2")


def bar_grad(a=None, b=None, c=None, d=None):
    """Used in multiple tests to simplify formatting of expected result"""
    ret = [("width", "10em")]
    if all(x is None for x in [a, b, c, d]):
        return ret
    return ret + [
        (
            "background",
            f"linear-gradient(90deg,{','.join([x for x in [a, b, c, d] if x])})",
        )
    ]


def no_bar():
    return bar_grad()


def bar_to(x, color="#d65f5f"):
    return bar_grad(f" {color} {x:.1f}%", f" transparent {x:.1f}%")


def bar_from_to(x, y, color="#d65f5f"):
    return bar_grad(
        f" transparent {x:.1f}%",
        f" {color} {x:.1f}%",
        f" {color} {y:.1f}%",
        f" transparent {y:.1f}%",
    )


class TestStylerBarAlign:
    def test_bar_align_mid_vmin(self):
        df = DataFrame({"A": [0, 1], "B": [-2, 4]})
        result = df.style.bar(align="mid", axis=None, vmin=-6)._compute().ctx
        expected = {
            (0, 0): bar_grad(),
            (1, 0): bar_grad(
                " transparent 60.0%",
                " #d65f5f 60.0%",
                " #d65f5f 70.0%",
                " transparent 70.0%",
            ),
            (0, 1): bar_grad(
                " transparent 40.0%",
                " #d65f5f 40.0%",
                " #d65f5f 60.0%",
                " transparent 60.0%",
            ),
            (1, 1): bar_grad(
                " transparent 60.0%",
                " #d65f5f 60.0%",
                " #d65f5f 100.0%",
                " transparent 100.0%",
            ),
        }
        assert result == expected

    def test_bar_align_mid_vmax(self):
        df = DataFrame({"A": [0, 1], "B": [-2, 4]})
        result = df.style.bar(align="mid", axis=None, vmax=8)._compute().ctx
        expected = {
            (0, 0): bar_grad(),
            (1, 0): bar_grad(
                " transparent 20.0%",
                " #d65f5f 20.0%",
                " #d65f5f 30.0%",
                " transparent 30.0%",
            ),
            (0, 1): bar_grad(
                " #d65f5f 20.0%",
                " transparent 20.0%",
            ),
            (1, 1): bar_grad(
                " transparent 20.0%",
                " #d65f5f 20.0%",
                " #d65f5f 60.0%",
                " transparent 60.0%",
            ),
        }
        assert result == expected

    def test_bar_align_mid_vmin_vmax_wide(self):
        df = DataFrame({"A": [0, 1], "B": [-2, 4]})
        result = df.style.bar(align="mid", axis=None, vmin=-3, vmax=7)._compute().ctx
        expected = {
            (0, 0): bar_grad(),
            (1, 0): bar_grad(
                " transparent 30.0%",
                " #d65f5f 30.0%",
                " #d65f5f 40.0%",
                " transparent 40.0%",
            ),
            (0, 1): bar_grad(
                " transparent 10.0%",
                " #d65f5f 10.0%",
                " #d65f5f 30.0%",
                " transparent 30.0%",
            ),
            (1, 1): bar_grad(
                " transparent 30.0%",
                " #d65f5f 30.0%",
                " #d65f5f 70.0%",
                " transparent 70.0%",
            ),
        }
        assert result == expected

    def test_bar_align_mid_vmin_vmax_clipping(self):
        df = DataFrame({"A": [0, 1], "B": [-2, 4]})
        result = df.style.bar(align="mid", axis=None, vmin=-1, vmax=3)._compute().ctx
        expected = {
            (0, 0): bar_grad(),
            (1, 0): bar_grad(
                " transparent 25.0%",
                " #d65f5f 25.0%",
                " #d65f5f 50.0%",
                " transparent 50.0%",
            ),
            (0, 1): bar_grad(" #d65f5f 25.0%", " transparent 25.0%"),
            (1, 1): bar_grad(
                " transparent 25.0%",
                " #d65f5f 25.0%",
                " #d65f5f 100.0%",
                " transparent 100.0%",
            ),
        }
        assert result == expected

    def test_bar_align_mid_nans(self):
        df = DataFrame({"A": [1, None], "B": [-1, 3]})
        result = df.style.bar(align="mid", axis=None)._compute().ctx
        expected = {
            (0, 0): bar_grad(
                " transparent 25.0%",
                " #d65f5f 25.0%",
                " #d65f5f 50.0%",
                " transparent 50.0%",
            ),
            (0, 1): bar_grad(" #d65f5f 25.0%", " transparent 25.0%"),
            (1, 0): bar_grad(),
            (1, 1): bar_grad(
                " transparent 25.0%",
                " #d65f5f 25.0%",
                " #d65f5f 100.0%",
                " transparent 100.0%",
            ),
        }
        assert result == expected

    def test_bar_align_zero_nans(self):
        df = DataFrame({"A": [1, None], "B": [-1, 2]})
        result = df.style.bar(align="zero", axis=None)._compute().ctx
        expected = {
            (0, 0): bar_grad(
                " transparent 50.0%",
                " #d65f5f 50.0%",
                " #d65f5f 75.0%",
                " transparent 75.0%",
            ),
            (0, 1): bar_grad(
                " transparent 25.0%",
                " #d65f5f 25.0%",
                " #d65f5f 50.0%",
                " transparent 50.0%",
            ),
            (1, 0): bar_grad(),
            (1, 1): bar_grad(
                " transparent 50.0%",
                " #d65f5f 50.0%",
                " #d65f5f 100.0%",
                " transparent 100.0%",
            ),
        }
        assert result == expected

    def test_bar_bad_align_raises(self):
        df = DataFrame({"A": [-100, -60, -30, -20]})
        msg = "`align` should be in {'left', 'right', 'mid', 'mean', 'zero'} or"
        with pytest.raises(ValueError, match=msg):
            df.style.bar(align="poorly", color=["#d65f5f", "#5fba7d"]).render()


@pytest.mark.parametrize(
    "align, exp",
    [
        ("left", [no_bar(), bar_to(50), bar_to(100)]),
        ("right", [bar_to(100), bar_from_to(50, 100), no_bar()]),
        ("mid", [bar_to(33.33), bar_to(66.66), bar_to(100)]),
        ("zero", [bar_from_to(50, 66.7), bar_from_to(50, 83.3), bar_from_to(50, 100)]),
        ("mean", [bar_to(50), no_bar(), bar_from_to(50, 100)]),
        (2.0, [bar_to(50), no_bar(), bar_from_to(50, 100)]),
        (np.median, [bar_to(50), no_bar(), bar_from_to(50, 100)]),
    ],
)
def test_align_positive_cases(align, exp):
    # test different align cases for all positive values
    data = DataFrame([[1], [2], [3]])
    result = data.style.bar(align=align)._compute().ctx
    expected = {(0, 0): exp[0], (1, 0): exp[1], (2, 0): exp[2]}
    assert result == expected


@pytest.mark.parametrize(
    "align, exp",
    [
        ("left", [bar_to(100), bar_to(50), no_bar()]),
        ("right", [no_bar(), bar_from_to(50, 100), bar_to(100)]),
        ("mid", [bar_from_to(66.66, 100), bar_from_to(33.33, 100), bar_to(100)]),
        ("zero", [bar_from_to(33.33, 50), bar_from_to(16.66, 50), bar_to(50)]),
        ("mean", [bar_from_to(50, 100), no_bar(), bar_to(50)]),
        (-2.0, [bar_from_to(50, 100), no_bar(), bar_to(50)]),
        (np.median, [bar_from_to(50, 100), no_bar(), bar_to(50)]),
    ],
)
def test_align_negative_cases(align, exp):
    # test different align cases for all negative values
    data = DataFrame([[-1], [-2], [-3]])
    result = data.style.bar(align=align)._compute().ctx
    expected = {(0, 0): exp[0], (1, 0): exp[1], (2, 0): exp[2]}
    assert result == expected


@pytest.mark.parametrize(
    "align, exp",
    [
        ("left", [no_bar(), bar_to(80), bar_to(100)]),
        ("right", [bar_to(100), bar_from_to(80, 100), no_bar()]),
        ("mid", [bar_to(60), bar_from_to(60, 80), bar_from_to(60, 100)]),
        ("zero", [bar_to(50), bar_from_to(50, 66.66), bar_from_to(50, 83.33)]),
        ("mean", [bar_to(50), bar_from_to(50, 66.66), bar_from_to(50, 83.33)]),
        (-0.0, [bar_to(50), bar_from_to(50, 66.66), bar_from_to(50, 83.33)]),
        (np.median, [bar_to(50), no_bar(), bar_from_to(50, 62.5)]),
    ],
)
def test_align_mixed_cases(align, exp):
    # test different align cases for mixed positive and negative values
    data = DataFrame([[-3], [1], [2]])
    result = data.style.bar(align=align)._compute().ctx
    expected = {(0, 0): exp[0], (1, 0): exp[1], (2, 0): exp[2]}
    assert result == expected


@pytest.mark.parametrize(
    "align, exp",
    [
        (
            "left",
            {
                "index": [[no_bar(), no_bar()], [bar_to(100), bar_to(100)]],
                "columns": [[no_bar(), bar_to(100)], [no_bar(), bar_to(100)]],
                "none": [[no_bar(), bar_to(33.33)], [bar_to(66.66), bar_to(100)]],
            },
        ),
        (
            "mid",
            {
                "index": [[bar_to(33.33), bar_to(50)], [bar_to(100), bar_to(100)]],
                "columns": [[bar_to(50), bar_to(100)], [bar_to(75), bar_to(100)]],
                "none": [[bar_to(25), bar_to(50)], [bar_to(75), bar_to(100)]],
            },
        ),
        (
            "zero",
            {
                "index": [
                    [bar_from_to(50, 66.66), bar_from_to(50, 75)],
                    [bar_from_to(50, 100), bar_from_to(50, 100)],
                ],
                "columns": [
                    [bar_from_to(50, 75), bar_from_to(50, 100)],
                    [bar_from_to(50, 87.5), bar_from_to(50, 100)],
                ],
                "none": [
                    [bar_from_to(50, 62.5), bar_from_to(50, 75)],
                    [bar_from_to(50, 87.5), bar_from_to(50, 100)],
                ],
            },
        ),
        (
            2,
            {
                "index": [
                    [bar_to(50), no_bar()],
                    [bar_from_to(50, 100), bar_from_to(50, 100)],
                ],
                "columns": [
                    [bar_to(50), no_bar()],
                    [bar_from_to(50, 75), bar_from_to(50, 100)],
                ],
                "none": [
                    [bar_from_to(25, 50), no_bar()],
                    [bar_from_to(50, 75), bar_from_to(50, 100)],
                ],
            },
        ),
    ],
)
@pytest.mark.parametrize("axis", ["index", "columns", "none"])
def test_align_axis(align, exp, axis):
    # test all axis combinations with positive values and different aligns
    data = DataFrame([[1, 2], [3, 4]])
    result = (
        data.style.bar(align=align, axis=None if axis == "none" else axis)
        ._compute()
        .ctx
    )
    expected = {
        (0, 0): exp[axis][0][0],
        (0, 1): exp[axis][0][1],
        (1, 0): exp[axis][1][0],
        (1, 1): exp[axis][1][1],
    }
    assert result == expected


def test_numerics():
    # test data is pre-selected for numeric values
    data = DataFrame([[1, "a"], [2, "b"]])
    result = data.style.bar()._compute().ctx
    assert (0, 1) not in result
    assert (1, 1) not in result


@pytest.mark.parametrize(
    "align, exp",
    [
        ("left", [no_bar(), bar_to(100, "green")]),
        ("right", [bar_to(100, "red"), no_bar()]),
        ("mid", [bar_to(25, "red"), bar_from_to(25, 100, "green")]),
        ("zero", [bar_from_to(33.33, 50, "red"), bar_from_to(50, 100, "green")]),
    ],
)
def test_colors_mixed(align, exp):
    data = DataFrame([[-1], [3]])
    result = data.style.bar(align=align, color=["red", "green"])._compute().ctx
    assert result == {(0, 0): exp[0], (1, 0): exp[1]}
