import pytest

from pandas import Index, MultiIndex


class TestIndexConstructor:
    # Tests for the Index constructor, specifically for cases that do
    #  not return a subclass

    def test_constructor_corner(self):
        # corner case
        msg = (
            r"Index\(\.\.\.\) must be called with a collection of some "
            "kind, 0 was passed"
        )
        with pytest.raises(TypeError, match=msg):
            Index(0)

    @pytest.mark.parametrize("index_vals", [[("A", 1), "B"], ["B", ("A", 1)]])
    def test_construction_list_mixed_tuples(self, index_vals):
        # see gh-10697: if we are constructing from a mixed list of tuples,
        # make sure that we are independent of the sorting order.
        index = Index(index_vals)
        assert isinstance(index, Index)
        assert not isinstance(index, MultiIndex)

    def test_constructor_wrong_kwargs(self):
        # GH #19348
        with pytest.raises(TypeError, match="Unexpected keyword arguments {'foo'}"):
            Index([], foo="bar")

    @pytest.mark.xfail(reason="see GH#21311: Index doesn't enforce dtype argument")
    def test_constructor_cast(self):
        msg = "could not convert string to float"
        with pytest.raises(ValueError, match=msg):
            Index(["a", "b", "c"], dtype=float)
