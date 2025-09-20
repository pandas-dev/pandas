# Copyright 2017 Google Inc.
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

"""Utilities for Google Media Downloads and Resumable Uploads.

===========
Subpackages
===========

Each subpackage is tailored to a specific transport library:

* the :mod:`~google.cloud.storage._media.requests` subpackage uses the ``requests``
  transport library.

.. _requests: http://docs.python-requests.org/
"""

from google.cloud.storage._media.common import UPLOAD_CHUNK_SIZE


__all__ = [
    "UPLOAD_CHUNK_SIZE",
]
