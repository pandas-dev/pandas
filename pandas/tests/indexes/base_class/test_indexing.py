import pytest

from pandas import Index


class TestGetSliceBounds:
    @pytest.mark.parametrize("kind", ["getitem", "loc", None])
    @pytest.mark.parametrize("side, expected", [("left", 4), ("right", 5)])
    def test_get_slice_bounds_within(self, kind, side, expected):
        index = Index(list("abcdef"))
        result = index.get_slice_bound("e", kind=kind, side=side)
        assert result == expected

    @pytest.mark.parametrize("kind", ["getitem", "loc", None])
    @pytest.mark.parametrize("side", ["left", "right"])
    @pytest.mark.parametrize(
        "data, bound, expected", [(list("abcdef"), "x", 6), (list("bcdefg"), "a", 0)]
    )
    def test_get_slice_bounds_outside(self, kind, side, expected, data, bound):
        index = Index(data)
        result = index.get_slice_bound(bound, kind=kind, side=side)
        assert result == expected

    def test_get_slice_bounds_invalid_side(self):
        with pytest.raises(ValueError, match="Invalid value for side kwarg"):
            Index([]).get_slice_bound("a", kind=None, side="middle")
