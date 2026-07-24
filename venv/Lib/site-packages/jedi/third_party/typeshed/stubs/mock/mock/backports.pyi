import sys
from unittest import IsolatedAsyncioTestCase as IsolatedAsyncioTestCase

if sys.version_info >= (3, 10):
    from inspect import iscoroutinefunction as iscoroutinefunction
else:
    from asyncio import iscoroutinefunction as iscoroutinefunction
