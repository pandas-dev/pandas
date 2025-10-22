"""
Tests for iloc dictionary assignment bug fix.

Regression test for GH#62723: Series.iloc assignment with dtype="object"
incorrectly converts dictionary values to Series.
"""

import pytest

from pandas import (
    DataFrame,
    Series,
)


class TestIlocDictionaryAssignment:
    """Tests for iloc dictionary assignment bug fix (GH#62723)."""

    def test_iloc_preserves_dict_object_dtype(self):
        """Test that iloc preserves dictionaries in object dtype Series."""
        # GH#62723
        s = Series(0, dtype="object")
        test_dict = {}

        # Assign dictionary via iloc
        s.iloc[0] = test_dict

        # Verify dictionary is preserved
        assert s[0] == test_dict
        assert isinstance(s[0], dict)

    def test_iloc_preserves_complex_dict_object_dtype(self):
        """Test that iloc preserves complex dictionaries in object dtype Series."""
        s = Series(0, dtype="object")
        test_dict = {
            "a": 1,
            "b": [1, 2, 3],
            "c": {"nested": True},
            "d": None,
            "e": 3.14,
        }

        s.iloc[0] = test_dict

        assert s[0] == test_dict
        assert isinstance(s[0], dict)
        assert s[0]["a"] == 1
        assert s[0]["b"] == [1, 2, 3]
        assert s[0]["c"] == {"nested": True}
        assert s[0]["d"] is None
        assert s[0]["e"] == 3.14

    def test_iloc_vs_loc_dict_assignment_consistency(self):
        """Test that iloc and direct assignment behave consistently."""
        # Original bug: s[0] = {} works but s.iloc[0] = {} converts to Series
        s = Series(0, dtype="object")

        # Direct assignment should work (baseline)
        s[0] = {}
        assert s[0] == {}
        assert isinstance(s[0], dict)

        # Reset and test iloc
        s = Series(0, dtype="object")
        s.iloc[0] = {}
        assert s[0] == {}
        assert isinstance(s[0], dict)

    def test_iloc_multiple_dict_assignments(self):
        """Test iloc dictionary assignment to multiple positions."""
        s = Series([0, 1, 2], dtype="object")
        dict1 = {"first": 1}
        dict2 = {"second": 2}

        s.iloc[0] = dict1
        s.iloc[2] = dict2

        assert s.iloc[0] == dict1
        assert isinstance(s.iloc[0], dict)
        assert s.iloc[1] == 1  # unchanged
        assert s.iloc[2] == dict2
        assert isinstance(s.iloc[2], dict)

    def test_iloc_dict_assignment_non_object_dtype_fails(self):
        """Test that iloc dict assignment to non-object dtypes fails as expected."""
        s = Series([1, 2, 3], dtype="int64")

        # This should fail for non-object dtypes
        with pytest.raises((ValueError, TypeError)):
            s.iloc[0] = {"key": "value"}

    def test_iloc_preserves_other_object_types(self):
        """Test that iloc preserves other object types correctly."""
        s = Series([None] * 4, dtype="object")

        # Test various object types
        test_list = [1, 2, 3]
        test_tuple = (1, 2, 3)
        test_set = {1, 2, 3}
        test_str = "test string"

        s.iloc[0] = test_list
        s.iloc[1] = test_tuple
        s.iloc[2] = test_set
        s.iloc[3] = test_str

        assert s.iloc[0] == test_list
        assert isinstance(s.iloc[0], list)
        assert s.iloc[1] == test_tuple
        assert isinstance(s.iloc[1], tuple)
        assert s.iloc[2] == test_set
        assert isinstance(s.iloc[2], set)
        assert s.iloc[3] == test_str
        assert isinstance(s.iloc[3], str)

    def test_dataframe_iloc_dict_assignment_unaffected(self):
        """Test that the fix doesn't break DataFrame iloc dict assignment."""
        df = DataFrame({"A": [0, 1], "B": [2, 3]}, dtype="object")
        test_dict = {"frame": "test"}

        # DataFrame iloc should still work correctly
        df.iloc[0, 0] = test_dict

        assert df.iloc[0, 0] == test_dict
        assert isinstance(df.iloc[0, 0], dict)

    def test_nested_dict_assignment(self):
        """Test iloc assignment with deeply nested dictionaries."""
        s = Series([None], dtype="object")
        nested_dict = {
            "level1": {
                "level2": {
                    "level3": {
                        "deep": "value",
                        "list": [1, 2, {"nested_in_list": True}],
                    }
                }
            },
            "another_key": [{"dict_in_list": "test"}],
        }

        s.iloc[0] = nested_dict

        assert s.iloc[0] == nested_dict
        assert isinstance(s.iloc[0], dict)
        assert s.iloc[0]["level1"]["level2"]["level3"]["deep"] == "value"
        assert s.iloc[0]["another_key"][0]["dict_in_list"] == "test"

    def test_original_bug_reproduction(self):
        """Direct test of the original bug report scenario."""
        # This is the exact code from the bug report
        s = Series(0, dtype="object")

        s[0] = {}
        assert s[0] == {}  # This should pass

        s.iloc[0] = {}
        assert s[0] == {}  # This was failing before the fix

        # Additional verification
        assert isinstance(s[0], dict)
        assert type(s[0]) == dict
