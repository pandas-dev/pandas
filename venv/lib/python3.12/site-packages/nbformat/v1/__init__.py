"""The main module for the v1 notebook format."""

# -----------------------------------------------------------------------------
#  Copyright (C) 2008-2011  The IPython Development Team
#
#  Distributed under the terms of the BSD License.  The full license is in
#  the file LICENSE, distributed as part of this software.
# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
# Imports
# -----------------------------------------------------------------------------
from __future__ import annotations

from .convert import upgrade
from .nbbase import NotebookNode, new_code_cell, new_notebook, new_text_cell
from .nbjson import reads as read_json
from .nbjson import reads as reads_json
from .nbjson import to_notebook as to_notebook_json
from .nbjson import writes as write_json
from .nbjson import writes as writes_json
