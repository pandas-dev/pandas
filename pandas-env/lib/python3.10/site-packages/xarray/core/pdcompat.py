# For reference, here is a copy of the pandas copyright notice:

# (c) 2011-2012, Lambda Foundry, Inc. and PyData Development Team
# All rights reserved.

# Copyright (c) 2008-2011 AQR Capital Management, LLC
# All rights reserved.

# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are
# met:

#     * Redistributions of source code must retain the above copyright
#        notice, this list of conditions and the following disclaimer.

#     * Redistributions in binary form must reproduce the above
#        copyright notice, this list of conditions and the following
#        disclaimer in the documentation and/or other materials provided
#        with the distribution.

#     * Neither the name of the copyright holder nor the names of any
#        contributors may be used to endorse or promote products derived
#        from this software without specific prior written permission.

# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDER AND CONTRIBUTORS
# "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
# LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
# A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
# OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
# SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
# LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
# DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
# THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
# (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
# OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
from __future__ import annotations

from enum import Enum
from typing import Literal

import pandas as pd
from packaging.version import Version


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


def nanosecond_precision_timestamp(*args, **kwargs) -> pd.Timestamp:
    """Return a nanosecond-precision Timestamp object.

    Note this function should no longer be needed after addressing GitHub issue
    #7493.
    """
    if Version(pd.__version__) >= Version("2.0.0"):
        return pd.Timestamp(*args, **kwargs).as_unit("ns")
    else:
        return pd.Timestamp(*args, **kwargs)
