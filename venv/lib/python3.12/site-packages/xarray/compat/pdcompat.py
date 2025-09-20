# For reference, here is a copy of the pandas copyright notice:

# BSD 3-Clause License

# Copyright (c) 2008-2011, AQR Capital Management, LLC, Lambda Foundry, Inc. and PyData Development Team
# All rights reserved.

# Copyright (c) 2011-2025, Open source contributors.

# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:

# * Redistributions of source code must retain the above copyright notice, this
#   list of conditions and the following disclaimer.

# * Redistributions in binary form must reproduce the above copyright notice,
#   this list of conditions and the following disclaimer in the documentation
#   and/or other materials provided with the distribution.

# * Neither the name of the copyright holder nor the names of its
#   contributors may be used to endorse or promote products derived from
#   this software without specific prior written permission.

# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
# DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
# FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
# DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
# SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
# CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
# OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
# OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

from __future__ import annotations

from enum import Enum
from typing import Literal

import pandas as pd

from xarray.core.types import PDDatetimeUnitOptions


def count_not_none(*args) -> int:
    """Compute the number of non-None arguments.

    Copied from pandas.core.common.count_not_none (not part of the public API)
    """
    return sum(arg is not None for arg in args)


class _NoDefault(Enum):
    """Used by pandas to specify a default value for a deprecated argument.
    Copied from pandas._libs.lib._NoDefault.

    See also:
    - pandas-dev/pandas#30788
    - pandas-dev/pandas#40684
    - pandas-dev/pandas#40715
    - pandas-dev/pandas#47045
    """

    no_default = "NO_DEFAULT"

    def __repr__(self) -> str:
        return "<no_default>"


no_default = (
    _NoDefault.no_default
)  # Sentinel indicating the default value following pandas
NoDefault = Literal[_NoDefault.no_default]  # For typing following pandas


def timestamp_as_unit(date: pd.Timestamp, unit: PDDatetimeUnitOptions) -> pd.Timestamp:
    """Convert the underlying int64 representation to the given unit.

    Compatibility function for pandas issue where "as_unit" is not defined
    for pandas.Timestamp in pandas versions < 2.2. Can be removed minimum
    pandas version is >= 2.2.
    """
    if hasattr(date, "as_unit"):
        date = date.as_unit(unit)
    elif hasattr(date, "_as_unit"):
        date = date._as_unit(unit)
    return date


def default_precision_timestamp(*args, **kwargs) -> pd.Timestamp:
    """Return a Timestamp object with the default precision.

    Xarray default is "ns".
    """
    dt = pd.Timestamp(*args, **kwargs)
    if dt.unit != "ns":
        dt = timestamp_as_unit(dt, "ns")
    return dt
