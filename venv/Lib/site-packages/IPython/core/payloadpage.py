"""A payload based version of page."""

# Copyright (c) IPython Development Team.
# Distributed under the terms of the Modified BSD License.

import warnings
from IPython.core.getipython import get_ipython

# see https://github.com/ipython/ipykernel/issues/1304
# this should be moved to ipykernel and removed in the long run.


def page(strng, start=0, screen_lines=0, pager_cmd=None):
    """Print a string, piping through a pager.

    This version ignores the screen_lines and pager_cmd arguments and uses
    IPython's payload system instead.

    Parameters
    ----------
    strng : str or mime-dict
        Text to page, or a mime-type keyed dict of already formatted data.
    start : int
        Starting line at which to place the display.
    """

    # Some routines may auto-compute start offsets incorrectly and pass a
    # negative value.  Offset to 0 for robustness.
    start = max(0, start)
    shell = get_ipython()

    if isinstance(strng, dict):
        data = strng
    else:
        data = {"text/plain": strng}
    payload = dict(
        source="page",
        data=data,
        start=start,
    )
    shell.payload_manager.write_payload(payload)
