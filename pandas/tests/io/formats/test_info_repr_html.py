import pytest

import pandas as pd


class Testfloatformat:
    @pytest.mark.parametrize(
        "data, format_option, expected_values",
        [
            ({"A": [12345.6789]}, "{:12.3f}", "12345.679"),
            ({"A": [None]}, "{:.3f}", "None"),
            ({"A": [""]}, "{:.2f}", ""),
            ({"A": [112345.6789]}, "{:6.3f}", "112345.679"),
        ],
    )  # test cases
    def test_float_formatting_html_output(self, data, format_option, expected_values):
        # set float format, avoid for string checks
        if format_option is not None:
            pd.set_option("display.float_format", format_option.format)

        # create dataframe
        df = pd.DataFrame(data)

        # capture html output
        html_output = df._repr_html_()

        # check
        assert expected_values in html_output

        # reset option
        if format_option is not None:
            pd.reset_option("display.float_format")
