#  Licensed to the Apache Software Foundation (ASF) under one
#  or more contributor license agreements.  See the NOTICE file
#  distributed with this work for additional information
#  regarding copyright ownership.  The ASF licenses this file
#  to you under the Apache License, Version 2.0 (the
#  "License"); you may not use this file except in compliance
#  with the License.  You may obtain a copy of the License at
#
#    http://www.apache.org/licenses/LICENSE-2.0
#
#  Unless required by applicable law or agreed to in writing,
#  software distributed under the License is distributed on an
#  "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY
#  KIND, either express or implied.  See the License for the
#  specific language governing permissions and limitations
#  under the License.
import functools
import warnings
from typing import Any, Callable, Optional


def deprecated(deprecated_in: str, removed_in: str, help_message: Optional[str] = None) -> Callable:  # type: ignore
    """Mark functions as deprecated.

    Adding this will result in a warning being emitted when the function is used.
    """
    if help_message is not None:
        help_message = f" {help_message}."

    def decorator(func: Callable):  # type: ignore
        @functools.wraps(func)
        def new_func(*args: Any, **kwargs: Any) -> Any:
            message = f"Call to {func.__name__}, deprecated in {deprecated_in}, will be removed in {removed_in}.{help_message}"

            _deprecation_warning(message)

            return func(*args, **kwargs)

        return new_func

    return decorator


def deprecation_notice(deprecated_in: str, removed_in: str, help_message: Optional[str]) -> str:
    """Return a deprecation notice."""
    return f"Deprecated in {deprecated_in}, will be removed in {removed_in}. {help_message}"


def deprecation_message(deprecated_in: str, removed_in: str, help_message: Optional[str]) -> None:
    """Mark properties or behaviors as deprecated.

    Adding this will result in a warning being emitted.
    """
    _deprecation_warning(deprecation_notice(deprecated_in, removed_in, help_message))


def _deprecation_warning(message: str) -> None:
    with warnings.catch_warnings():  # temporarily override warning handling
        warnings.warn(
            message,
            category=DeprecationWarning,
            stacklevel=2,
        )
