from io import StringIO

import pytest


@pytest.mark.parametrize(
    "kwds",
    [{"on_bad_lines": "warn"}, {"error_bad_lines": True}, {"warn_bad_lines": True}],
)
def test_bad_lines_fails(pyarrow_parser_only, kwds):
    # GH#
    data = """a,b,c
1,2,3
"""
    parser = pyarrow_parser_only
    msg = (
        f"The '{list(kwds.keys())[0]}' option is not supported "
        f"with the 'pyarrow' engine"
    )
    with pytest.raises(ValueError, match=msg):
        parser.read_csv(StringIO(data), **kwds)
