# Copyright 2016 Google LLC
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#      http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

"""Google Auth Library for Python."""

import logging
import sys
import warnings

from google.auth import version as google_auth_version
from google.auth._default import (
    default,
    load_credentials_from_dict,
    load_credentials_from_file,
)


__version__ = google_auth_version.__version__


__all__ = ["default", "load_credentials_from_file", "load_credentials_from_dict"]


class Python37DeprecationWarning(DeprecationWarning):  # pragma: NO COVER
    """
    Deprecation warning raised when Python 3.7 runtime is detected.
    Python 3.7 support will be dropped after January 1, 2024.
    """

    pass


# Raise warnings for deprecated versions
eol_message = (
    "You are using a Python version {} past its end of life. Google will update "
    "google-auth with critical bug fixes on a best-effort basis, but not "
    "with any other fixes or features. Please upgrade your Python version, "
    "and then update google-auth."
)
if sys.version_info.major == 3 and sys.version_info.minor == 8:  # pragma: NO COVER
    warnings.warn(eol_message.format("3.8"), FutureWarning)
elif sys.version_info.major == 3 and sys.version_info.minor == 9:  # pragma: NO COVER
    warnings.warn(eol_message.format("3.9"), FutureWarning)

# Set default logging handler to avoid "No handler found" warnings.
logging.getLogger(__name__).addHandler(logging.NullHandler())
