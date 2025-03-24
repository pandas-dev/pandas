# Licensed to the Apache Software Foundation (ASF) under one
# or more contributor license agreements.  See the NOTICE file
# distributed with this work for additional information
# regarding copyright ownership.  The ASF licenses this file
# to you under the Apache License, Version 2.0 (the
# "License"); you may not use this file except in compliance
# with the License.  You may obtain a copy of the License at
#
#   http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing,
# software distributed under the License is distributed on an
# "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY
# KIND, either express or implied.  See the License for the
# specific language governing permissions and limitations
# under the License.

# cython: language_level = 3

"""
For debugging, install crash handlers that print a backtrace.
"""

import threading

cdef extern from "backward.hpp" nogil:
    cdef struct CSignalHandling"backward::SignalHandling":
        pass


cdef class _SignalHandling:
    cdef CSignalHandling _c_signal_handler


_CRASH_HANDLER = None
_CRASH_HANDLER_LOCK = threading.Lock()


def _install_crash_handler():
    global _CRASH_HANDLER
    with _CRASH_HANDLER_LOCK:
        if _CRASH_HANDLER:
            return
        _CRASH_HANDLER = _SignalHandling()
