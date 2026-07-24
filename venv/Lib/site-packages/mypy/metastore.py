"""Interfaces for accessing metadata.

We provide two implementations.
 * The "classic" file system implementation, which uses a directory
   structure of files.
 * A hokey sqlite backed implementation, which basically simulates
   the file system in an effort to work around poor file system performance
   on OS X.
"""

from __future__ import annotations

import binascii
import os
import time
from abc import abstractmethod
from collections.abc import Iterable
from typing import TYPE_CHECKING, Any

from mypy.util import hash_path_stem, os_path_join

if TYPE_CHECKING:
    # We avoid importing sqlite3 unless we are using it so we can mostly work
    # on semi-broken pythons that are missing it.
    import sqlite3


class MetadataStore:
    """Generic interface for metadata storage."""

    @abstractmethod
    def getmtime(self, name: str) -> float:
        """Read the mtime of a metadata entry.

        Raises FileNotFound if the entry does not exist.
        """

    @abstractmethod
    def read(self, name: str) -> bytes:
        """Read the contents of a metadata entry.

        Raises FileNotFound if the entry does not exist.
        """

    @abstractmethod
    def write(self, name: str, data: bytes, mtime: float | None = None) -> bool:
        """Write a metadata entry.

        If mtime is specified, set it as the mtime of the entry. Otherwise,
        the current time is used.

        Returns True if the entry is successfully written, False otherwise.
        """

    @abstractmethod
    def remove(self, name: str) -> None:
        """Delete a metadata entry"""

    @abstractmethod
    def commit(self) -> None:
        """If the backing store requires a commit, do it.

        But N.B. that this is not *guaranteed* to do anything, and
        there is no guarantee that changes are not made until it is
        called.
        """

    def commit_path(self, name: str) -> None:
        """Commit changes related to a specific cache path.

        For sharded stores, this commits only the shard containing the path.
        Default implementation commits everything.
        """
        self.commit()

    @abstractmethod
    def list_all(self) -> Iterable[str]: ...

    @abstractmethod
    def close(self) -> None:
        """Release any resources held by the backing store."""


def random_string() -> str:
    return binascii.hexlify(os.urandom(8)).decode("ascii")


class FilesystemMetadataStore(MetadataStore):
    def __init__(self, cache_dir_prefix: str) -> None:
        # We check startswith instead of equality because the version
        # will have already been appended by the time the cache dir is
        # passed here.
        if cache_dir_prefix.startswith(os.devnull):
            self.cache_dir_prefix = None
        else:
            self.cache_dir_prefix = cache_dir_prefix

    def getmtime(self, name: str) -> float:
        if not self.cache_dir_prefix:
            raise FileNotFoundError()

        return int(os.path.getmtime(os_path_join(self.cache_dir_prefix, name)))

    def read(self, name: str) -> bytes:
        assert not os.path.isabs(name), "Don't use absolute paths!"

        if not self.cache_dir_prefix:
            raise FileNotFoundError()

        with open(os_path_join(self.cache_dir_prefix, name), "rb", buffering=0) as f:
            return f.read()

    def write(self, name: str, data: bytes, mtime: float | None = None) -> bool:
        assert not os.path.isabs(name), "Don't use absolute paths!"

        if not self.cache_dir_prefix:
            return False

        path = os_path_join(self.cache_dir_prefix, name)
        tmp_filename = path + "." + random_string()
        try:
            os.makedirs(os.path.dirname(path), exist_ok=True)
            with open(tmp_filename, "wb") as f:
                f.write(data)
            os.replace(tmp_filename, path)
            if mtime is not None:
                os.utime(path, times=(mtime, mtime))

        except OSError:
            return False
        return True

    def remove(self, name: str) -> None:
        if not self.cache_dir_prefix:
            raise FileNotFoundError()

        os.remove(os_path_join(self.cache_dir_prefix, name))

    def commit(self) -> None:
        pass

    def list_all(self) -> Iterable[str]:
        if not self.cache_dir_prefix:
            return

        for dir, _, files in os.walk(self.cache_dir_prefix):
            dir = os.path.relpath(dir, self.cache_dir_prefix)
            for file in files:
                yield os.path.normpath(os_path_join(dir, file))

    def close(self) -> None:
        pass


SCHEMA = """
CREATE TABLE IF NOT EXISTS files2 (
    path TEXT UNIQUE NOT NULL,
    mtime REAL,
    data BLOB
);
CREATE INDEX IF NOT EXISTS path_idx on files2(path);
"""


def connect_db(db_file: str, set_journal_mode: bool) -> sqlite3.Connection:
    import sqlite3.dbapi2

    db = sqlite3.dbapi2.connect(db_file, check_same_thread=False)
    # This is a bit unfortunate (as we may get corrupt cache after e.g. Ctrl + C),
    # but without this flag, commits are *very* slow, especially when using HDDs,
    # see https://www.sqlite.org/faq.html#q19 for details.
    db.execute("PRAGMA synchronous=OFF")
    if set_journal_mode:
        db.execute("PRAGMA journal_mode=WAL")
    db.executescript(SCHEMA)
    return db


class SqliteMetadataStore(MetadataStore):
    def __init__(
        self, cache_dir_prefix: str, set_journal_mode: bool = False, num_shards: int = 1
    ) -> None:
        # We check startswith instead of equality because the version
        # will have already been appended by the time the cache dir is
        # passed here.
        self.dbs: list[sqlite3.Connection] = []
        self.num_shards = num_shards
        self.dirty_shards: set[int] = set()
        if cache_dir_prefix.startswith(os.devnull):
            return

        os.makedirs(cache_dir_prefix, exist_ok=True)
        if num_shards <= 1:
            self.dbs.append(
                connect_db(os_path_join(cache_dir_prefix, "cache.db"), set_journal_mode)
            )
        else:
            for i in range(num_shards):
                self.dbs.append(
                    connect_db(os_path_join(cache_dir_prefix, f"cache.{i}.db"), set_journal_mode)
                )

    def _shard_index(self, name: str) -> int:
        if self.num_shards <= 1:
            return 0
        return hash_path_stem(name) % self.num_shards

    def _db_for(self, name: str) -> sqlite3.Connection:
        if not self.dbs:
            raise FileNotFoundError()
        return self.dbs[self._shard_index(name)]

    def _query(self, name: str, field: str) -> Any:
        # Raises FileNotFound for consistency with the file system version
        db = self._db_for(name)
        cur = db.execute(f"SELECT {field} FROM files2 WHERE path = ?", (name,))
        results = cur.fetchall()
        if not results:
            raise FileNotFoundError()
        assert len(results) == 1
        return results[0][0]

    def getmtime(self, name: str) -> float:
        mtime = self._query(name, "mtime")
        assert isinstance(mtime, float)
        return mtime

    def read(self, name: str) -> bytes:
        data = self._query(name, "data")
        assert isinstance(data, bytes)
        return data

    def write(self, name: str, data: bytes, mtime: float | None = None) -> bool:
        import sqlite3

        if not self.dbs:
            return False
        try:
            if mtime is None:
                mtime = time.time()
            db = self._db_for(name)
            db.execute(
                "INSERT OR REPLACE INTO files2(path, mtime, data) VALUES(?, ?, ?)",
                (name, mtime, data),
            )
            self.dirty_shards.add(self._shard_index(name))
        except sqlite3.OperationalError:
            return False
        return True

    def remove(self, name: str) -> None:
        db = self._db_for(name)
        db.execute("DELETE FROM files2 WHERE path = ?", (name,))
        self.dirty_shards.add(self._shard_index(name))

    def commit(self) -> None:
        for i in self.dirty_shards:
            self.dbs[i].commit()
        self.dirty_shards.clear()

    def commit_path(self, name: str) -> None:
        i = self._shard_index(name)
        if i in self.dirty_shards:
            self.dbs[i].commit()
            self.dirty_shards.discard(i)

    def list_all(self) -> Iterable[str]:
        for db in self.dbs:
            for row in db.execute("SELECT path FROM files2"):
                yield row[0]

    def close(self) -> None:
        for db in self.dbs:
            db.close()
        self.dbs.clear()

    def __del__(self) -> None:
        self.close()
