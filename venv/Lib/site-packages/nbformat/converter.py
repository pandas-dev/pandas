"""API for converting notebooks between versions."""

# Copyright (c) IPython Development Team.
# Distributed under the terms of the Modified BSD License.
from __future__ import annotations

from . import versions
from .reader import get_version
from .validator import ValidationError


def convert(nb, to_version):
    """Convert a notebook node object to a specific version.  Assumes that
    all the versions starting from 1 to the latest major X are implemented.
    In other words, there should never be a case where v1 v2 v3 v5 exist without
    a v4.  Also assumes that all conversions can be made in one step increments
    between major versions and ignores minor revisions.

    Parameters
    ----------
    nb : NotebookNode
    to_version : int
        Major revision to convert the notebook to.  Can either be an upgrade or
        a downgrade.

    Raises
    ------
    ValueError
        Notebook failed to convert.
    ValueError
        The version specified is invalid or doesn't exist.
    ValidationError
        Conversion failed due to missing expected attributes.
    """

    # Get input notebook version.
    (version, version_minor) = get_version(nb)

    # Check if destination is target version, if so return contents
    if version == to_version:
        return nb

    # If the version exist, try to convert to it one step at a time.
    if to_version in versions:
        # Get the the version that this recursion will convert to as a step
        # closer to the final revision.  Make sure the newer of the conversion
        # functions is used to perform the conversion.
        if to_version > version:
            step_version = version + 1
            convert_function = versions[step_version].upgrade
        else:
            step_version = version - 1
            convert_function = versions[version].downgrade

        try:
            # Convert and make sure version changed during conversion.
            converted = convert_function(nb)
            if converted.get("nbformat", 1) == version:
                msg = "Failed to convert notebook from v%d to v%d." % (version, step_version)
                raise ValueError(msg)
        except AttributeError as e:
            msg = f"Notebook could not be converted from version {version} to version {step_version} because it's missing a key: {e}"
            raise ValidationError(msg) from None

        # Recursively convert until target version is reached.
        return convert(converted, to_version)
    raise ValueError(
        "Cannot convert notebook to v%d because that version doesn't exist" % (to_version)
    )
