# This file is part of Hypothesis, which may be found at
# https://github.com/HypothesisWorks/hypothesis/
#
# Copyright the Hypothesis Authors.
# Individual contributors are listed in AUTHORS.rst and the git log.
#
# This Source Code Form is subject to the terms of the Mozilla Public License,
# v. 2.0. If a copy of the MPL was not distributed with this file, You can
# obtain one at https://mozilla.org/MPL/2.0/.

import datetime
import warnings

from hypothesis.errors import HypothesisDeprecationWarning


def note_deprecation(
    message: str, *, since: str, has_codemod: bool, stacklevel: int = 0
) -> None:
    if since != "RELEASEDAY":
        date = datetime.date.fromisoformat(since)
        assert datetime.date(2021, 1, 1) <= date
    if has_codemod:
        message += (
            "\n    The `hypothesis codemod` command-line tool can automatically "
            "refactor your code to fix this warning."
        )
    warnings.warn(HypothesisDeprecationWarning(message), stacklevel=2 + stacklevel)
