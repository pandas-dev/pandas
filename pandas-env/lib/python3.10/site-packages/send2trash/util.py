# encoding: utf-8
# Copyright 2017 Virgil Dupras

# This software is licensed under the "BSD" License as described in the "LICENSE" file,
# which should be included with this package. The terms are also available at
# http://www.hardcoded.net/licenses/bsd_license

from send2trash.compat import text_type, binary_type, iterable_type


def preprocess_paths(paths):
    if isinstance(paths, iterable_type) and not isinstance(paths, (text_type, binary_type)):
        paths = list(paths)
    elif not isinstance(paths, list):
        paths = [paths]
    # Convert items such as pathlib paths to strings
    paths = [path.__fspath__() if hasattr(path, "__fspath__") else path for path in paths]
    return paths
