"""
Tests for metadata preservation in groupby operations.
"""

import numpy as np

import pandas._testing as tm


class TestGroupByMetadataPreservation:
    def test_groupby_apply_preserves_metadata(self):
        """Test that groupby.apply() preserves _metadata from subclassed DataFrame."""
        # Create a subclassed DataFrame with metadata
        subdf = tm.SubclassedDataFrame(
            {"X": [1, 1, 2, 2, 3], "Y": np.arange(0, 5), "Z": np.arange(10, 15)}
        )
        subdf.testattr = "test"

        # Apply groupby operation
        result = subdf.groupby("X").apply(np.sum, axis=0, include_groups=False)

        # Check that metadata is preserved
        assert hasattr(result, "testattr"), (
            "Metadata attribute 'testattr' should be preserved"
        )
        assert result.testattr == "test", "Metadata value should be preserved"

        # Compare with equivalent operation that preserves metadata
        expected = subdf.groupby("X").sum()
        assert expected.testattr == "test", (
            "Equivalent operation should preserve metadata"
        )
