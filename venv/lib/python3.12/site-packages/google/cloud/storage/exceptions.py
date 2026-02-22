# Copyright 2024 Google LLC
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

"""Exceptions raised by the library."""

# These exceptions were originally part of the google-resumable-media library
# but were integrated into python-storage in version 3.0. For backwards
# compatibility with applications which use except blocks with
# google-resumable-media exceptions, if the library google-resumable-media is
# installed, make all exceptions subclasses of the exceptions from that library.
# Note that either way, the classes will subclass Exception, either directly or
# indirectly.
#
# This backwards compatibility feature may be removed in a future major version
# update. Please update application code to use the new exception classes in
# this module.
try:
    from google.resumable_media import InvalidResponse as InvalidResponseDynamicParent
    from google.resumable_media import DataCorruption as DataCorruptionDynamicParent
except ImportError:
    InvalidResponseDynamicParent = Exception
    DataCorruptionDynamicParent = Exception


class InvalidResponse(InvalidResponseDynamicParent):
    """Error class for responses which are not in the correct state.

    Args:
        response (object): The HTTP response which caused the failure.
        args (tuple): The positional arguments typically passed to an
            exception class.
    """

    def __init__(self, response, *args):
        if InvalidResponseDynamicParent is Exception:
            super().__init__(*args)
            self.response = response
            """object: The HTTP response object that caused the failure."""
        else:
            super().__init__(response, *args)


class DataCorruption(DataCorruptionDynamicParent):
    """Error class for corrupt media transfers.

    Args:
        response (object): The HTTP response which caused the failure.
        args (tuple): The positional arguments typically passed to an
            exception class.
    """

    def __init__(self, response, *args):
        if DataCorruptionDynamicParent is Exception:
            super().__init__(*args)
            self.response = response
            """object: The HTTP response object that caused the failure."""
        else:
            super().__init__(response, *args)
