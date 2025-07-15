# For reference, here is a copy of the scikit-learn copyright notice:

# BSD 3-Clause License

# Copyright (c) 2007-2021 The scikit-learn developers.
# All rights reserved.

# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:

# * Redistributions of source code must retain the above copyright notice, this
#   list of conditions and the following disclaimer.

# * Redistributions in binary form must reproduce the above copyright notice,
#   this list of conditions and the following disclaimer in the documentation
#   and/or other materials provided with the distribution.

# * Neither the name of the copyright holder nor the names of its
#   contributors may be used to endorse or promote products derived from
#   this software without specific prior written permission.

# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
# DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
# FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
# DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
# SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
# CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
# OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
# OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE


import inspect
import warnings
from collections.abc import Callable
from functools import wraps
from typing import TypeVar

from xarray.core.utils import emit_user_level_warning

T = TypeVar("T", bound=Callable)

POSITIONAL_OR_KEYWORD = inspect.Parameter.POSITIONAL_OR_KEYWORD
KEYWORD_ONLY = inspect.Parameter.KEYWORD_ONLY
POSITIONAL_ONLY = inspect.Parameter.POSITIONAL_ONLY
EMPTY = inspect.Parameter.empty


def _deprecate_positional_args(version) -> Callable[[T], T]:
    """Decorator for methods that issues warnings for positional arguments

    Using the keyword-only argument syntax in pep 3102, arguments after the
    ``*`` will issue a warning when passed as a positional argument.

    Parameters
    ----------
    version : str
        version of the library when the positional arguments were deprecated

    Examples
    --------
    Deprecate passing `b` as positional argument:

    def func(a, b=1):
        pass

    @_deprecate_positional_args("v0.1.0")
    def func(a, *, b=2):
        pass

    func(1, 2)

    Notes
    -----
    This function is adapted from scikit-learn under the terms of its license. See
    licences/SCIKIT_LEARN_LICENSE
    """

    def _decorator(func):
        signature = inspect.signature(func)

        pos_or_kw_args = []
        kwonly_args = []
        for name, param in signature.parameters.items():
            if param.kind in (POSITIONAL_OR_KEYWORD, POSITIONAL_ONLY):
                pos_or_kw_args.append(name)
            elif param.kind == KEYWORD_ONLY:
                kwonly_args.append(name)
                if param.default is EMPTY:
                    # IMHO `def f(a, *, b):` does not make sense -> disallow it
                    # if removing this constraint -> need to add these to kwargs as well
                    raise TypeError("Keyword-only param without default disallowed.")

        @wraps(func)
        def inner(*args, **kwargs):
            name = func.__name__
            n_extra_args = len(args) - len(pos_or_kw_args)
            if n_extra_args > 0:
                extra_args = ", ".join(kwonly_args[:n_extra_args])

                warnings.warn(
                    f"Passing '{extra_args}' as positional argument(s) to {name} "
                    f"was deprecated in version {version} and will raise an error two "
                    "releases later. Please pass them as keyword arguments."
                    "",
                    FutureWarning,
                    stacklevel=2,
                )

                zip_args = zip(
                    kwonly_args[:n_extra_args], args[-n_extra_args:], strict=True
                )
                kwargs.update(zip_args)

                return func(*args[:-n_extra_args], **kwargs)

            return func(*args, **kwargs)

        return inner

    return _decorator


def deprecate_dims(func: T, old_name="dims") -> T:
    """
    For functions that previously took `dims` as a kwarg, and have now transitioned to
    `dim`. This decorator will issue a warning if `dims` is passed while forwarding it
    to `dim`.
    """

    @wraps(func)
    def wrapper(*args, **kwargs):
        if old_name in kwargs:
            emit_user_level_warning(
                f"The `{old_name}` argument has been renamed to `dim`, and will be removed "
                "in the future. This renaming is taking place throughout xarray over the "
                "next few releases.",
                # Upgrade to `DeprecationWarning` in the future, when the renaming is complete.
                PendingDeprecationWarning,
            )
            kwargs["dim"] = kwargs.pop(old_name)
        return func(*args, **kwargs)

    # We're quite confident we're just returning `T` from this function, so it's fine to ignore typing
    # within the function.
    return wrapper  # type: ignore[return-value]
