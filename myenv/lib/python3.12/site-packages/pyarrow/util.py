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

# Miscellaneous utility code

import os
import contextlib
import functools
import gc
import socket
import sys
import textwrap
import types
import warnings


_DEPR_MSG = (
    "pyarrow.{} is deprecated as of {}, please use pyarrow.{} instead."
)


def doc(*docstrings, **params):
    """
    A decorator that takes docstring templates, concatenates them, and finally
    performs string substitution on them.
    This decorator will add a variable "_docstring_components" to the wrapped
    callable to keep track of the original docstring template for potential future use.
    If the docstring is a template, it will be saved as a string.
    Otherwise, it will be saved as a callable and the docstring will be obtained via
    the __doc__ attribute.
    This decorator cannot be used on Cython classes due to a CPython constraint,
    which enforces the __doc__ attribute to be read-only.
    See https://github.com/python/cpython/issues/91309

    Parameters
    ----------
    *docstrings : None, str, or callable
        The string / docstring / docstring template to be prepended in order
        before the default docstring under the callable.
    **params
        The key/value pairs used to format the docstring template.
    """

    def decorator(decorated):
        docstring_components = []

        # collect docstrings and docstring templates
        for docstring in docstrings:
            if docstring is None:
                continue
            if hasattr(docstring, "_docstring_components"):
                docstring_components.extend(
                    docstring._docstring_components
                )
            elif isinstance(docstring, str) or docstring.__doc__:
                docstring_components.append(docstring)

        # append the callable's docstring last
        if decorated.__doc__:
            docstring_components.append(textwrap.dedent(decorated.__doc__))

        params_applied = [
            component.format(**params)
            if isinstance(component, str) and len(params) > 0
            else component
            for component in docstring_components
        ]

        decorated.__doc__ = "".join(
            [
                component
                if isinstance(component, str)
                else textwrap.dedent(component.__doc__ or "")
                for component in params_applied
            ]
        )

        decorated._docstring_components = (
            docstring_components
        )
        return decorated

    return decorator


def _deprecate_api(old_name, new_name, api, next_version, type=FutureWarning):
    msg = _DEPR_MSG.format(old_name, next_version, new_name)

    def wrapper(*args, **kwargs):
        warnings.warn(msg, type)
        return api(*args, **kwargs)
    return wrapper


def _deprecate_class(old_name, new_class, next_version,
                     instancecheck=True):
    """
    Raise warning if a deprecated class is used in an isinstance check.
    """
    class _DeprecatedMeta(type):
        def __instancecheck__(self, other):
            warnings.warn(
                _DEPR_MSG.format(old_name, next_version, new_class.__name__),
                FutureWarning,
                stacklevel=2
            )
            return isinstance(other, new_class)

    return _DeprecatedMeta(old_name, (new_class,), {})


def _is_iterable(obj):
    try:
        iter(obj)
        return True
    except TypeError:
        return False


def _is_path_like(path):
    return isinstance(path, str) or hasattr(path, '__fspath__')


def _stringify_path(path):
    """
    Convert *path* to a string or unicode path if possible.
    """
    if isinstance(path, str):
        return os.path.expanduser(path)

    # checking whether path implements the filesystem protocol
    try:
        return os.path.expanduser(path.__fspath__())
    except AttributeError:
        pass

    raise TypeError("not a path-like object")


def product(seq):
    """
    Return a product of sequence items.
    """
    return functools.reduce(lambda a, b: a*b, seq, 1)


def get_contiguous_span(shape, strides, itemsize):
    """
    Return a contiguous span of N-D array data.

    Parameters
    ----------
    shape : tuple
    strides : tuple
    itemsize : int
      Specify array shape data

    Returns
    -------
    start, end : int
      The span end points.
    """
    if not strides:
        start = 0
        end = itemsize * product(shape)
    else:
        start = 0
        end = itemsize
        for i, dim in enumerate(shape):
            if dim == 0:
                start = end = 0
                break
            stride = strides[i]
            if stride > 0:
                end += stride * (dim - 1)
            elif stride < 0:
                start += stride * (dim - 1)
        if end - start != itemsize * product(shape):
            raise ValueError('array data is non-contiguous')
    return start, end


def find_free_port():
    sock = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
    with contextlib.closing(sock) as sock:
        sock.bind(('', 0))
        sock.setsockopt(socket.SOL_SOCKET, socket.SO_REUSEADDR, 1)
        return sock.getsockname()[1]


def guid():
    from uuid import uuid4
    return uuid4().hex


def _break_traceback_cycle_from_frame(frame):
    # Clear local variables in all inner frames, so as to break the
    # reference cycle.
    this_frame = sys._getframe(0)
    refs = gc.get_referrers(frame)
    while refs:
        for frame in refs:
            if frame is not this_frame and isinstance(frame, types.FrameType):
                break
        else:
            # No frame found in referrers (finished?)
            break
        refs = None
        # Clear the frame locals, to try and break the cycle (it is
        # somewhere along the chain of execution frames).
        frame.clear()
        # To visit the inner frame, we need to find it among the
        # referrers of this frame (while `frame.f_back` would let
        # us visit the outer frame).
        refs = gc.get_referrers(frame)
    refs = frame = this_frame = None


def download_tzdata_on_windows():
    r"""
    Download and extract latest IANA timezone database into the
    location expected by Arrow which is %USERPROFILE%\Downloads\tzdata.
    """
    if sys.platform != 'win32':
        raise TypeError(f"Timezone database is already provided by {sys.platform}")

    import tarfile

    tzdata_path = os.path.expandvars(r"%USERPROFILE%\Downloads\tzdata")
    tzdata_compressed = os.path.join(tzdata_path, "tzdata.tar.gz")
    os.makedirs(tzdata_path, exist_ok=True)

    from urllib.request import urlopen
    with urlopen('https://data.iana.org/time-zones/tzdata-latest.tar.gz') as response:
        with open(tzdata_compressed, 'wb') as f:
            f.write(response.read())

    assert os.path.exists(tzdata_compressed)

    tarfile.open(tzdata_compressed).extractall(tzdata_path)

    with urlopen('https://raw.githubusercontent.com/unicode-org/cldr/master/common/supplemental/windowsZones.xml') as response_zones:   # noqa
        with open(os.path.join(tzdata_path, "windowsZones.xml"), 'wb') as f:
            f.write(response_zones.read())
