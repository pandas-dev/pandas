import numpy as np
import pytest

from pandas.util._validators import validate_percentile


def test_validate_percentile_negative_array_input():
    percentile_array = np.array([-0.5, 0.25, 1.5])
    expected_error_msg = "percentiles should all be in the interval [0, 1]"
    with pytest.raises(ValueError, match=(expected_error_msg)):
        validate_percentile(percentile_array)
