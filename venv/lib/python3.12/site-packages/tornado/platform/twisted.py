# Licensed under the Apache License, Version 2.0 (the "License"); you may
# not use this file except in compliance with the License. You may obtain
# a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS, WITHOUT
# WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. See the
# License for the specific language governing permissions and limitations
# under the License.
"""Bridges between the Twisted package and Tornado."""

import sys

from twisted.internet.defer import Deferred  # type: ignore
from twisted.python import failure  # type: ignore

from tornado.concurrent import Future, future_set_exc_info
from tornado import gen

import typing  # noqa: F401


def install() -> None:
    """Install ``AsyncioSelectorReactor`` as the default Twisted reactor.

    .. deprecated:: 5.1

       This function is provided for backwards compatibility; code
       that does not require compatibility with older versions of
       Tornado should use
       ``twisted.internet.asyncioreactor.install()`` directly.

    .. versionchanged:: 6.0.3

       In Tornado 5.x and before, this function installed a reactor
       based on the Tornado ``IOLoop``. When that reactor
       implementation was removed in Tornado 6.0.0, this function was
       removed as well. It was restored in Tornado 6.0.3 using the
       ``asyncio`` reactor instead.

    """
    from twisted.internet.asyncioreactor import install  # type: ignore

    install()


if hasattr(gen.convert_yielded, "register"):

    @gen.convert_yielded.register(Deferred)
    def _(d: Deferred) -> Future:
        f = Future()  # type: Future[typing.Any]

        def errback(failure: failure.Failure) -> None:
            try:
                failure.raiseException()
                # Should never happen, but just in case
                raise Exception("errback called without error")
            except:
                future_set_exc_info(f, sys.exc_info())

        d.addCallbacks(f.set_result, errback)
        return f
