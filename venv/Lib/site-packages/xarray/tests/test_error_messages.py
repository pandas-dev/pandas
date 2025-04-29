"""
This new file is intended to test the quality & friendliness of error messages that are
raised by xarray. It's currently separate from the standard tests, which are more
focused on the functions working (though we could consider integrating them.).
"""

import pytest


def test_no_var_in_dataset(ds):
    with pytest.raises(
        KeyError,
        match=(
            r"No variable named 'foo'. Variables on the dataset include \['z1', 'z2', 'x', 'time', 'c', 'y'\]"
        ),
    ):
        ds["foo"]
