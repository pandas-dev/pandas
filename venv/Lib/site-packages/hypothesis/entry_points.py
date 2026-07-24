# This file is part of Hypothesis, which may be found at
# https://github.com/HypothesisWorks/hypothesis/
#
# Copyright the Hypothesis Authors.
# Individual contributors are listed in AUTHORS.rst and the git log.
#
# This Source Code Form is subject to the terms of the Mozilla Public License,
# v. 2.0. If a copy of the MPL was not distributed with this file, You can
# obtain one at https://mozilla.org/MPL/2.0/.

"""Run all functions registered for the "hypothesis" entry point.

This can be used with `st.register_type_strategy` to register strategies for your
custom types, running the relevant code when *hypothesis* is imported instead of
your package.
"""

import importlib.metadata
import os


def run() -> None:
    if os.environ.get("HYPOTHESIS_NO_PLUGINS"):
        return

    for entry in importlib.metadata.entry_points(
        group="hypothesis"
    ):  # pragma: no cover
        hook = entry.load()
        if callable(hook):
            hook()
