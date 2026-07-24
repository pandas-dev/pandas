"""Translation related utilities. When imported, injects _ to builtins"""

# Copyright (c) Jupyter Development Team.
# Distributed under the terms of the Modified BSD License.
import gettext
import os
import warnings


def _trans_gettext_deprecation_helper(*args, **kwargs):
    """The trans gettext deprecation helper."""
    warn_msg = "The alias `_()` will be deprecated. Use `_i18n()` instead."
    warnings.warn(warn_msg, FutureWarning, stacklevel=2)
    return trans.gettext(*args, **kwargs)


# Set up message catalog access
base_dir = os.path.realpath(os.path.join(__file__, "..", ".."))
trans = gettext.translation(
    "notebook", localedir=os.path.join(base_dir, "notebook/i18n"), fallback=True
)
_ = _trans_gettext_deprecation_helper
_i18n = trans.gettext
