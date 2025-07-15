# Author: Ovidiu Predescu
# Date: July 2011
#
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

import sys
import unittest

from tornado.testing import AsyncTestCase, gen_test

try:
    from twisted.internet.defer import inlineCallbacks  # type: ignore

    have_twisted = True
except ImportError:
    have_twisted = False
except Exception:
    # Twisted is currently incompatible with the first 3.14 alpha release; disable this
    # test until the beta when it will hopefully be fixed (note that this requires us to
    # update our requirements.txt to pick up a new version of twisted).
    if sys.version_info[:2] == (3, 14) and sys.version_info.releaselevel == "alpha":
        have_twisted = False
    else:
        raise
else:
    # Not used directly but needed for `yield deferred` to work.
    import tornado.platform.twisted  # noqa: F401

skipIfNoTwisted = unittest.skipUnless(have_twisted, "twisted module not present")


@skipIfNoTwisted
class ConvertDeferredTest(AsyncTestCase):
    @gen_test
    def test_success(self):
        @inlineCallbacks
        def fn():
            if False:
                # inlineCallbacks doesn't work with regular functions;
                # must have a yield even if it's unreachable.
                yield
            return 42

        res = yield fn()
        self.assertEqual(res, 42)

    @gen_test
    def test_failure(self):
        @inlineCallbacks
        def fn():
            if False:
                yield
            1 / 0

        with self.assertRaises(ZeroDivisionError):
            yield fn()


if __name__ == "__main__":
    unittest.main()
