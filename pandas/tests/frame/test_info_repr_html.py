import pytest
import pandas as pd


class Testfloatformat:
    @pytest.mark.parametrize("data, format_option, expected_values", [
        ({"A": [12345.6789]}, "{:12.3f}", "&nbsp;&nbsp;&nbsp;12345.679"),
        ({"A": [None]}, "{:.3f}", "&nbsp;None"),
        ({"A": [""]}, "{:.2f}", "&nbsp;"),
        ({"A": [112345.6789]}, "{:6.3f}", "112345.679"),
        ({"A": ["foo    foo"]}, None, "&nbsp;foo&nbsp;&nbsp;&nbsp;&nbsp;foo"),
        ({"A": [None]}, None, "&nbsp;None"),
        ({"A": ["foo foo foo"]}, None, "&nbsp;foo&nbsp;foo&nbsp;foo"),
    ]) #test cases

    def test_float_formatting_html_output(self,data,format_option, expected_values):
            # set float format, avoid for string checks
            if format_option is not None:
                 pd.set_option("display.float_format", format_option.format)
            
            # crate dataframe
            df = pd.DataFrame(data)
                              
            # capture html output
            html_output = df._repr_html_()

            # check
            assert expected_values in html_output

            # reset option
            if format_option is not None:
                 pd.reset_option("display.float_format")