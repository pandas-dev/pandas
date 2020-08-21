import pandas as pd


class TestFlags:
    def test_equality(self):
        a = pd.DataFrame().set_flags(allows_duplicate_labels=True).flags
        b = pd.DataFrame().set_flags(allows_duplicate_labels=False).flags

        assert a == a
        assert b == b
        assert a != b
        assert a != 2

    def test_set(self):
        df = pd.DataFrame().set_flags(allows_duplicate_labels=True)
        a = df.flags
        a.allows_duplicate_labels = False
        assert a.allows_duplicate_labels is False
        a["allows_duplicate_labels"] = True
        assert a.allows_duplicate_labels is True

    def test_repr(self):
        a = repr(pd.DataFrame({"A"}).set_flags(allows_duplicate_labels=True).flags)
        assert a == "<Flags(allows_duplicate_labels=True)>"
        a = repr(pd.DataFrame({"A"}).set_flags(allows_duplicate_labels=False).flags)
        assert a == "<Flags(allows_duplicate_labels=False)>"
