"""Directory-listing cache (``dircache``) update strategies.

The filesystem keeps a cache of directory listings in ``self.dircache``. After
any mutating operation (delete, write, move, ...) that cache must be kept
consistent with the bucket. This module groups that concern into two mixins so
the bookkeeping lives in one place rather than being interleaved with the I/O
in :mod:`gcsfs.core` and :mod:`gcsfs.extended_gcsfs`:

- :class:`DirCacheUpdater` -- the default strategy. For deletes and moves it
  broadly invalidates the affected parent directory and all of its ancestors, so
  the next listing re-reads them from the bucket. For writes it invalidates only
  the immediate parent when that parent is already cached (and therefore known to
  exist), falling back to broad ancestor invalidation otherwise.
- :class:`HnsDirCacheUpdater` -- a targeted strategy for Hierarchical Namespace
  (HNS) buckets, where directories are first-class persistent objects. It
  mutates parent listings in place for deletes and moves, falling back to the
  broad strategy for non-HNS paths.

Both are mixed into a filesystem class and rely on the host class for
``dircache``, ``split_path``, ``_parent``, ``_strip_protocol``,
``invalidate_cache``, ``_is_bucket_hns_enabled`` and ``_process_object``.
"""

import asyncio


class DirCacheUpdater:
    """Default ``dircache`` update strategy.

    Deletes and moves broadly invalidate the affected parent and all of its
    ancestors; writes invalidate only the immediate parent when it is already
    cached (and therefore known to exist), otherwise they fall back to broad
    ancestor invalidation.

    Mixed into :class:`gcsfs.core.GCSFileSystem`. The methods are the hooks that
    the mutating I/O paths (``_rm_file``, ``_rm_files``, ``_put_file``, ...) call
    after a successful operation; subclasses override them to provide a more
    targeted strategy.
    """

    async def _mv_file_cache_update(self, path1, path2, response=None):
        self.invalidate_cache(self._parent(path1))
        self.invalidate_cache(self._parent(path2))

    async def _rm_file_cache_update(self, path):
        # Single-path form of _rm_files_cache_update; kept as one code path so
        # subclasses only need to override the batch method.
        await self._rm_files_cache_update([path])

    async def _rm_files_cache_update(self, paths):
        parents = set(self._parent(p) for p in paths) | set(paths)
        for parent in parents:
            self.invalidate_cache(parent)

    async def _write_file_cache_update(self, path):
        # A file was created or overwritten at ``path`` (via put/pipe/cp).
        #
        # When the immediate parent's listing is already cached, we have listed
        # it before, so it already exists and adding a file inside it cannot
        # create a new directory in any ancestor's listing. Invalidate only that
        # immediate parent so the new file is picked up on its next listing.
        #
        # When the parent is not cached, this write may have implicitly created
        # it (and intermediate directories) -- flat buckets simulate directories
        # from object prefixes, and HNS buckets auto-create missing parents on
        # object write. Leaving a cached ancestor untouched would hide the new
        # directory, so fall back to invalidating the parent and all ancestors.
        #
        # Like the rest of this module, this relies on the dircache being
        # consistent with the bucket (this client is the sole mutator between
        # listings); concurrent external mutations are reconciled only by
        # ``invalidate_cache`` / listings expiry, not here.
        parent = self._parent(path)
        if parent in self.dircache:
            self.dircache.pop(parent, None)
        else:
            self.invalidate_cache(parent)


class HnsDirCacheUpdater(DirCacheUpdater):
    """Targeted ``dircache`` update strategy for HNS buckets.

    Mixed into :class:`gcsfs.extended_gcsfs.ExtendedGcsFileSystem`. For HNS
    buckets directories are persistent objects, so deletes and moves usually
    only change the contents of the immediate parent's listing. Non-HNS /
    cross-bucket paths defer to :class:`DirCacheUpdater` via ``super()``. The
    write path is bucket-type-agnostic and handled entirely by
    :meth:`DirCacheUpdater._write_file_cache_update`.
    """

    def _cache_drop_entries(self, parent, names):
        """Remove entries whose ``name`` is in ``names`` from ``parent``'s
        cached listing. No-op when ``parent`` is not cached.

        ``names`` is a set of stripped ``bucket/key`` names matching the
        ``name`` field of the cached entries.
        """
        if parent in self.dircache:
            self.dircache[parent] = [
                e for e in self.dircache[parent] if e.get("name") not in names
            ]

    def _cache_add_entry(self, parent, entry):
        """Append ``entry`` to ``parent``'s cached listing. No-op when
        ``parent`` is not cached."""
        if parent in self.dircache:
            self.dircache[parent].append(entry)

    def _cache_upsert_entry(self, parent, entry):
        """Replace any cached entry with the same name before appending
        ``entry``. No-op when ``parent`` is not cached."""
        if parent in self.dircache:
            name = entry.get("name")
            self.dircache[parent] = [
                e for e in self.dircache[parent] if e.get("name") != name
            ]
            self.dircache[parent].append(entry)

    @staticmethod
    def _directory_cache_entry(name, key):
        """Build a cached directory-listing entry for an HNS folder located at
        ``name`` (stripped ``bucket/key`` path) with object key ``key``."""
        return {
            "Key": key,
            "Size": 0,
            "name": name,
            "size": 0,
            "type": "directory",
            "storageClass": "DIRECTORY",
        }

    def _update_dircache_after_rename(self, path1, path2):
        """
        Performs a targeted update of the directory cache after a successful
        folder rename operation.

        This involves three main steps:
        1. Removing the source folder and all its descendants from the cache.
        2. Removing the source folder's entry from its parent's listing.
        3. Adding the new destination folder's entry to its parent's listing.

        Args:
            path1 (str): The source path that was renamed.
            path2 (str): The destination path.
        """
        # dircache keys and entry names are stored without the protocol, so
        # normalize the incoming paths first; otherwise the pop/startswith
        # matching below silently misses every cached entry for a ``gs://`` path.
        path1 = self._strip_protocol(path1)
        path2 = self._strip_protocol(path2)

        # 1. Find and remove all descendant paths of the source from the cache.
        source_prefix = f"{path1.rstrip('/')}/"
        for key in list(self.dircache):
            if key.startswith(source_prefix):
                self.dircache.pop(key, None)

        # 2. Remove the old source entry from its parent's listing.
        self.dircache.pop(path1, None)
        self._cache_drop_entries(self._parent(path1), {path1})

        # 3. Invalidate the destination path/subtree and update its parent's cache.
        dest_prefix = f"{path2.rstrip('/')}/"
        for key in list(self.dircache):
            if key == path2 or key.startswith(dest_prefix):
                self.dircache.pop(key, None)
        _, key2, _ = self.split_path(path2)
        self._cache_upsert_entry(
            self._parent(path2), self._directory_cache_entry(path2, key2)
        )

    async def _mv_file_cache_update(self, path1, path2, response=None):
        """
        Update the cache after a file move operation.

        For HNS-enabled buckets where the move is within the same bucket, this method
        directly updates the directory cache by removing the source entry from it's
        parent cache and adding destination path as a new entry in it's corresponding parent cache.
        This avoids invalidating the entire parent directory cache, which is beneficial for HNS
        performance.

        For non-HNS buckets or cross-bucket moves, it falls back to the default
        behavior (invalidating the cache for both source and destination parents).
        """
        src_bucket, _, _ = self.split_path(path1)
        dest_bucket, _, _ = self.split_path(path2)

        if await self._is_bucket_hns_enabled(src_bucket) and src_bucket == dest_bucket:
            # Source: removing the entry never changes an ancestor's listing, so
            # drop it in place.
            self._cache_drop_entries(self._parent(path1), {self._strip_protocol(path1)})
            dest_parent = self._parent(path2)
            if response and dest_parent in self.dircache:
                # Destination parent is already cached (so it existed before the
                # move) -> updating it in place cannot leave an ancestor stale.
                self._cache_upsert_entry(
                    dest_parent, self._process_object(dest_bucket, response)
                )
            else:
                # The destination parent may have been created by this move;
                # broad-invalidate it so a cached ancestor can't hide the new
                # directory (also covers a missing/empty move response).
                self.invalidate_cache(dest_parent)
        else:
            await super()._mv_file_cache_update(path1, path2, response)

    async def _rm_files_cache_update(self, paths):
        """
        Update the cache after a (batch) file delete operation.

        For HNS-enabled buckets, directories are first-class persistent objects,
        so deleting a file only changes the contents of its immediate parent.
        Each deleted entry is removed from its immediate parent's listing
        (grouped per parent so each listing is rewritten only once), avoiding the
        redundant invalidation of all ancestor directories up to the root.

        Entries are matched by name only; the object generation is not taken
        into account. Non-HNS paths fall back to the default broad invalidation
        of the parents and all of their ancestors.
        """
        # Split each path once and resolve each distinct bucket's HNS status
        # once, concurrently, rather than awaiting a (cached) lookup per path.
        split_paths = [(path, *self.split_path(path)) for path in paths]
        buckets = list({bucket for _, bucket, _, _ in split_paths})
        hns_enabled = dict(
            zip(
                buckets,
                await asyncio.gather(
                    *(self._is_bucket_hns_enabled(bucket) for bucket in buckets)
                ),
            )
        )

        # Group the names to drop by their immediate parent so each listing is
        # rewritten only once. Non-HNS paths fall back to broad invalidation.
        removed_names_by_parent = {}
        non_hns_paths = []
        for path, bucket, key, _ in split_paths:
            if hns_enabled[bucket]:
                removed_names_by_parent.setdefault(self._parent(path), set()).add(
                    f"{bucket}/{key}"
                )
            else:
                non_hns_paths.append(path)

        for parent, removed_names in removed_names_by_parent.items():
            self._cache_drop_entries(parent, removed_names)

        if non_hns_paths:
            await super()._rm_files_cache_update(non_hns_paths)
