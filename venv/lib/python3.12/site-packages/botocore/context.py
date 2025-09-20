# Copyright 2025 Amazon.com, Inc. or its affiliates. All Rights Reserved.
#
# Licensed under the Apache License, Version 2.0 (the "License"). You
# may not use this file except in compliance with the License. A copy of
# the License is located at
#
# http://aws.amazon.com/apache2.0/
#
# or in the "license" file accompanying this file. This file is
# distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF
# ANY KIND, either express or implied. See the License for the specific
# language governing permissions and limitations under the License.
"""
NOTE: All classes and functions in this module are considered private and are
subject to abrupt breaking changes. Please do not use them directly.
"""

from contextlib import contextmanager
from contextvars import ContextVar
from copy import deepcopy
from dataclasses import dataclass, field
from functools import wraps


@dataclass
class ClientContext:
    """
    Encapsulation of objects tracked within the ``_context`` context variable.

    ``features`` is a set responsible for storing features used during
    preparation of an AWS request. ``botocore.useragent.register_feature_id``
    is used to add to this set.
    """

    features: set[str] = field(default_factory=set)


_context = ContextVar("_context")


def get_context():
    """Get the current ``_context`` context variable if set, else None."""
    return _context.get(None)


def set_context(ctx):
    """Set the current ``_context`` context variable.

    :type ctx: ClientContext
    :param ctx: Client context object to set as the current context variable.

    :rtype: contextvars.Token
    :returns: Token object used to revert the context variable to what it was
        before the corresponding set.
    """
    token = _context.set(ctx)
    return token


def reset_context(token):
    """Reset the current ``_context`` context variable.

    :type token: contextvars.Token
    :param token: Token object to reset the context variable.
    """
    _context.reset(token)


@contextmanager
def start_as_current_context(ctx=None):
    """
    Context manager that copies the passed or current context object and sets
    it as the current context variable. If no context is found, a new
    ``ClientContext`` object is created. It mainly ensures the context variable
    is reset to the previous value once the executed code returns.

    Example usage:

        def my_feature():
            with start_as_current_context():
                register_feature_id('MY_FEATURE')
                pass

    :type ctx: ClientContext
    :param ctx: The client context object to set as the new context variable.
        If not provided, the current or a new context variable is used.
    """
    current = ctx or get_context()
    if current is None:
        new = ClientContext()
    else:
        new = deepcopy(current)
    token = set_context(new)
    try:
        yield
    finally:
        reset_context(token)


def with_current_context(hook=None):
    """
    Decorator that wraps ``start_as_current_context`` and optionally invokes a
    hook within the newly-set context. This is just syntactic sugar to avoid
    indenting existing code under the context manager.

    Example usage:

        @with_current_context(partial(register_feature_id, 'MY_FEATURE'))
        def my_feature():
            pass

    :type hook: callable
    :param hook: A callable that will be invoked within the scope of the
        ``start_as_current_context`` context manager.
    """

    def decorator(func):
        @wraps(func)
        def wrapper(*args, **kwargs):
            with start_as_current_context():
                if hook:
                    hook()
                return func(*args, **kwargs)

        return wrapper

    return decorator
