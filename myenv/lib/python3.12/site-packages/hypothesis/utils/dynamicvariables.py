# This file is part of Hypothesis, which may be found at
# https://github.com/HypothesisWorks/hypothesis/
#
# Copyright the Hypothesis Authors.
# Individual contributors are listed in AUTHORS.rst and the git log.
#
# This Source Code Form is subject to the terms of the Mozilla Public License,
# v. 2.0. If a copy of the MPL was not distributed with this file, You can
# obtain one at https://mozilla.org/MPL/2.0/.

import threading
from contextlib import contextmanager


class DynamicVariable:
    def __init__(self, default):
        self.default = default
        self.data = threading.local()

    @property
    def value(self):
        return getattr(self.data, "value", self.default)

    @value.setter
    def value(self, value):
        self.data.value = value

    @contextmanager
    def with_value(self, value):
        old_value = self.value
        try:
            self.data.value = value
            yield
        finally:
            self.data.value = old_value
