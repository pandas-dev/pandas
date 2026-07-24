# encoding: utf-8
# Copyright 2017 Virgil Dupras

# This software is licensed under the "BSD" License as described in the "LICENSE" file,
# which should be included with this package. The terms are also available at
# http://www.hardcoded.net/licenses/bsd_license

import os
import collections.abc


def preprocess_paths(paths):
    if isinstance(paths, collections.abc.Iterable) and not isinstance(paths, (str, bytes)):
        paths = list(paths)
    else:
        paths = [paths]
    # Convert items such as pathlib paths to strings
    return [os.fspath(path) for path in paths]
