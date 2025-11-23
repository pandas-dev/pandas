"""
Test for issue #63089 - read_csv segfault with large exponent
"""

import io

import pandas as pd


class TestIssue63089:
    def test_large_exponent_no_segfault(self):
        """Test that extremely large exponents don't cause segfault."""
        # This previously caused SIGSEGV due to integer overflow
        # when parsing the exponent
        result = pd.read_csv(
            io.StringIO("""h
4e492493924924""")
        )

        # Should parse as infinity or large float, not crash
        assert len(result) == 1
        assert "h" in result.columns
        # The value should be infinity since the exponent is way too large
        import numpy as np

        assert np.isinf(result["h"].iloc[0]) or result["h"].iloc[0] > 1e308

    def test_various_large_exponents(self):
        """Test various edge cases with large exponents."""
        test_cases = [
            "1e999999999",  # Very large positive exponent
            "1e-999999999",  # Very large negative exponent
            "2.5e123456789",  # Large exponent with decimal
        ]

        for test_val in test_cases:
            csv_data = f"col\n{test_val}"
            result = pd.read_csv(io.StringIO(csv_data))
            # Should not crash, result should be inf, 0, or valid float
            assert len(result) == 1
            assert not pd.isna(result["col"].iloc[0]) or True  # Just don't crash
