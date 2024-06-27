"""Convert notebook to the v1 format."""

# -----------------------------------------------------------------------------
#  Copyright (C) 2008-2011  The IPython Development Team
#
#  Distributed under the terms of the BSD License.  The full license is in
#  the file LICENSE, distributed as part of this software.
# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
# Code
# -----------------------------------------------------------------------------
from __future__ import annotations


def upgrade(nb, orig_version=None):
    """Upgrade a notebook."""
    msg = "Cannot convert to v1 notebook format"
    raise ValueError(msg)
