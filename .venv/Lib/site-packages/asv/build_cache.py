# Licensed under a 3-clause BSD style license - see LICENSE.rst

import os

from . import util


class BuildCache:
    """
    Build cache

    Data is cached in a directory tree::

        {self._path}/
            {self._path}/{commit_hash}/*
            {self._path}/{commit_hash}.timestamp

    If the timestamp file is missing, the subdirectory is ignored (and
    subject to cleanup).

    The cache cleanup retains the latest ``build_cache_size`` items
    that have a valid timestamp file.

    The timestamp files are created by ``self.finalize_cache_dir(commmit_hash)``,
    which also triggers a cache cleanup.

    The finalization should be called only after package is installed successfully,
    keeping in mind that ``build_cache_size`` may be 0.

    """

    def __init__(self, conf, root):
        self._root = root
        self._path = os.path.join(root, 'asv-build-cache')
        self._cache_size = getattr(conf, 'build_cache_size', 2)

    def _get_cache_dir(self, commit_hash):
        """
        Get the cache dir and timestamp file corresponding to a given commit hash.
        """
        path = os.path.join(self._path, commit_hash)
        stamp = path + ".timestamp"
        return path, stamp

    def _remove_cache_dir(self, commit_hash):
        path, stamp = self._get_cache_dir(commit_hash)
        if os.path.isdir(path):
            util.long_path_rmtree(path)
        if os.path.exists(stamp):
            os.unlink(stamp)

    def _get_cache_contents(self):
        """
        Return list of commit hash directories in the cache (containing
        wheels or not), sorted by decreasing timestamp
        """
        if not os.path.isdir(self._path):
            return []

        def sort_key(name):
            path, stamp = self._get_cache_dir(name)
            try:
                return os.stat(stamp).st_mtime
            except OSError:
                return 0

        names = os.listdir(self._path)
        names.sort(key=sort_key, reverse=True)
        return names

    def _cleanup_build_cache(self):
        # First remove items without timestamp
        if os.path.isdir(self._path):
            names = os.listdir(self._path)
            for name in names:
                path, stamp = self._get_cache_dir(name)
                if not os.path.exists(stamp):
                    self._remove_cache_dir(name)

        # Then remove old items
        names = self._get_cache_contents()
        for name in names[self._cache_size:]:
            self._remove_cache_dir(name)

    def get_cache_dir(self, commit_hash):
        path, stamp = self._get_cache_dir(commit_hash)
        if (os.path.isdir(path) and
                os.path.isfile(stamp) and
                os.listdir(path)):
            return path

        return None

    def create_cache_dir(self, commit_hash):
        self._remove_cache_dir(commit_hash)

        path, stamp = self._get_cache_dir(commit_hash)
        os.makedirs(path)
        return path

    def finalize_cache_dir(self, commit_hash):
        path, stamp = self._get_cache_dir(commit_hash)

        if os.path.isdir(path) and os.listdir(path):
            # Finalize
            with open(stamp, 'wb'):
                pass

        # Cleanup build cache
        self._cleanup_build_cache()
