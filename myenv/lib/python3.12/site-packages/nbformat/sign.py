"""Utilities for signing notebooks"""

# Copyright (c) IPython Development Team.
# Distributed under the terms of the Modified BSD License.
from __future__ import annotations

import hashlib
import os
import sys
import typing as t
from collections import OrderedDict
from contextlib import contextmanager
from datetime import datetime, timezone
from hmac import HMAC
from pathlib import Path

try:
    import sqlite3

    # Use adapters recommended by Python 3.12 stdlib docs.
    # https://docs.python.org/3.12/library/sqlite3.html#default-adapters-and-converters-deprecated
    def adapt_datetime_iso(val):
        """Adapt datetime.datetime to timezone-naive ISO 8601 date."""
        return val.isoformat()

    def convert_datetime(val):
        """Convert ISO 8601 datetime to datetime.datetime object."""
        return datetime.fromisoformat(val.decode())

    sqlite3.register_adapter(datetime, adapt_datetime_iso)
    sqlite3.register_converter("datetime", convert_datetime)
except ImportError:
    try:
        from pysqlite2 import dbapi2 as sqlite3  # type:ignore[no-redef]
    except ImportError:
        sqlite3 = None  # type:ignore[assignment]

from base64 import encodebytes

from jupyter_core.application import JupyterApp, base_flags
from traitlets import Any, Bool, Bytes, Callable, Enum, Instance, Integer, Unicode, default, observe
from traitlets.config import LoggingConfigurable, MultipleInstanceError

from . import NO_CONVERT, __version__, read, reads

algorithms_set = hashlib.algorithms_guaranteed
# The shake algorithms in are not compatible with hmac
# due to required length argument in digests
algorithms = [a for a in algorithms_set if not a.startswith("shake_")]


class SignatureStore:
    """Base class for a signature store."""

    def store_signature(self, digest, algorithm):
        """Implement in subclass to store a signature.

        Should not raise if the signature is already stored.
        """
        raise NotImplementedError

    def check_signature(self, digest, algorithm):
        """Implement in subclass to check if a signature is known.

        Return True for a known signature, False for unknown.
        """
        raise NotImplementedError

    def remove_signature(self, digest, algorithm):
        """Implement in subclass to delete a signature.

        Should not raise if the signature is not stored.
        """
        raise NotImplementedError

    def close(self):
        """Close any open connections this store may use.

        If the store maintains any open connections (e.g. to a database),
        they should be closed.
        """


class MemorySignatureStore(SignatureStore):
    """Non-persistent storage of signatures in memory."""

    cache_size = 65535

    def __init__(self):
        """Initialize a memory signature store."""
        # We really only want an ordered set, but the stdlib has OrderedDict,
        # and it's easy to use a dict as a set.
        self.data = OrderedDict()

    def store_signature(self, digest, algorithm):
        """Store a signature."""
        key = (digest, algorithm)
        # Pop it so it goes to the end when we reinsert it
        self.data.pop(key, None)
        self.data[key] = None

        self._maybe_cull()

    def _maybe_cull(self):
        """If more than cache_size signatures are stored, delete the oldest 25%"""
        if len(self.data) < self.cache_size:
            return

        for _ in range(len(self.data) // 4):
            self.data.popitem(last=False)

    def check_signature(self, digest, algorithm):
        """Check a signature."""
        key = (digest, algorithm)
        if key in self.data:
            # Move it to the end (.move_to_end() method is new in Py3)
            del self.data[key]
            self.data[key] = None
            return True
        return False

    def remove_signature(self, digest, algorithm):
        """Remove a signature."""
        self.data.pop((digest, algorithm), None)


class SQLiteSignatureStore(SignatureStore, LoggingConfigurable):
    """Store signatures in an SQLite database."""

    # 64k entries ~ 12MB
    cache_size = Integer(
        65535,
        help="""The number of notebook signatures to cache.
        When the number of signatures exceeds this value,
        the oldest 25% of signatures will be culled.
        """,
    ).tag(config=True)

    def __init__(self, db_file, **kwargs):
        """Initialize a sql signature store."""
        super().__init__(**kwargs)
        self.db_file = db_file
        self.db = self._connect_db(db_file)

    def close(self):
        """Close the db."""
        if self.db is not None:
            self.db.close()

    def _connect_db(self, db_file):
        kwargs: dict[str, t.Any] = {
            "detect_types": sqlite3.PARSE_DECLTYPES | sqlite3.PARSE_COLNAMES
        }
        db = None
        try:
            db = sqlite3.connect(db_file, **kwargs)
            self.init_db(db)
        except (sqlite3.DatabaseError, sqlite3.OperationalError):
            if db_file != ":memory:":
                old_db_location = db_file + ".bak"
                if db is not None:
                    db.close()
                self.log.warning(
                    (
                        "The signatures database cannot be opened; maybe it is corrupted or encrypted. "
                        "You may need to rerun your notebooks to ensure that they are trusted to run Javascript. "
                        "The old signatures database has been renamed to %s and a new one has been created."
                    ),
                    old_db_location,
                )
                try:
                    Path(db_file).rename(old_db_location)
                    db = sqlite3.connect(db_file, **kwargs)
                    self.init_db(db)
                except (sqlite3.DatabaseError, sqlite3.OperationalError, OSError):
                    if db is not None:
                        db.close()
                    self.log.warning(
                        "Failed committing signatures database to disk. "
                        "You may need to move the database file to a non-networked file system, "
                        "using config option `NotebookNotary.db_file`. "
                        "Using in-memory signatures database for the remainder of this session."
                    )
                    self.db_file = ":memory:"
                    db = sqlite3.connect(":memory:", **kwargs)
                    self.init_db(db)
            else:
                raise
        return db

    def init_db(self, db):
        """Initialize the db."""
        db.execute(
            """
            CREATE TABLE IF NOT EXISTS nbsignatures
            (
                id integer PRIMARY KEY AUTOINCREMENT,
                algorithm text,
                signature text,
                path text,
                last_seen timestamp
            )"""
        )
        db.execute(
            """
            CREATE INDEX IF NOT EXISTS algosig ON nbsignatures(algorithm, signature)
            """
        )
        db.commit()

    def store_signature(self, digest, algorithm):
        """Store a signature in the db."""
        if self.db is None:
            return
        if not self.check_signature(digest, algorithm):
            self.db.execute(
                """
                INSERT INTO nbsignatures (algorithm, signature, last_seen)
                VALUES (?, ?, ?)
                """,
                (algorithm, digest, datetime.now(tz=timezone.utc)),
            )
        else:
            self.db.execute(
                """UPDATE nbsignatures SET last_seen = ? WHERE
                algorithm = ? AND
                signature = ?;
                """,
                (datetime.now(tz=timezone.utc), algorithm, digest),
            )
        self.db.commit()

        # Check size and cull old entries if necessary
        (n,) = self.db.execute("SELECT Count(*) FROM nbsignatures").fetchone()
        if n > self.cache_size:
            self.cull_db()

    def check_signature(self, digest, algorithm):
        """Check a signature against the db."""
        if self.db is None:
            return False
        r = self.db.execute(
            """SELECT id FROM nbsignatures WHERE
            algorithm = ? AND
            signature = ?;
            """,
            (algorithm, digest),
        ).fetchone()
        if r is None:
            return False
        self.db.execute(
            """UPDATE nbsignatures SET last_seen = ? WHERE
            algorithm = ? AND
            signature = ?;
            """,
            (datetime.now(tz=timezone.utc), algorithm, digest),
        )
        self.db.commit()
        return True

    def remove_signature(self, digest, algorithm):
        """Remove a signature from the db."""
        self.db.execute(
            """DELETE FROM nbsignatures WHERE
                algorithm = ? AND
                signature = ?;
            """,
            (algorithm, digest),
        )

        self.db.commit()

    def cull_db(self):
        """Cull oldest 25% of the trusted signatures when the size limit is reached"""
        self.db.execute(
            """DELETE FROM nbsignatures WHERE id IN (
            SELECT id FROM nbsignatures ORDER BY last_seen DESC LIMIT -1 OFFSET ?
        );
        """,
            (max(int(0.75 * self.cache_size), 1),),
        )


def yield_everything(obj):
    """Yield every item in a container as bytes

    Allows any JSONable object to be passed to an HMAC digester
    without having to serialize the whole thing.
    """
    if isinstance(obj, dict):
        for key in sorted(obj):
            value = obj[key]
            assert isinstance(key, str)
            yield key.encode()
            yield from yield_everything(value)
    elif isinstance(obj, (list, tuple)):
        for element in obj:
            yield from yield_everything(element)
    elif isinstance(obj, str):
        yield obj.encode("utf8")
    else:
        yield str(obj).encode("utf8")


def yield_code_cells(nb):
    """Iterator that yields all cells in a notebook

    nbformat version independent
    """
    if nb.nbformat >= 4:
        for cell in nb["cells"]:
            if cell["cell_type"] == "code":
                yield cell
    elif nb.nbformat == 3:
        for ws in nb["worksheets"]:
            for cell in ws["cells"]:
                if cell["cell_type"] == "code":
                    yield cell


@contextmanager
def signature_removed(nb):
    """Context manager for operating on a notebook with its signature removed

    Used for excluding the previous signature when computing a notebook's signature.
    """
    save_signature = nb["metadata"].pop("signature", None)
    try:
        yield
    finally:
        if save_signature is not None:
            nb["metadata"]["signature"] = save_signature


class NotebookNotary(LoggingConfigurable):
    """A class for computing and verifying notebook signatures."""

    data_dir = Unicode(help="""The storage directory for notary secret and database.""").tag(
        config=True
    )

    @default("data_dir")
    def _data_dir_default(self):
        app = None
        try:
            if JupyterApp.initialized():
                app = JupyterApp.instance()
        except MultipleInstanceError:
            pass
        if app is None:
            # create an app, without the global instance
            app = JupyterApp()
            app.initialize(argv=[])
        return app.data_dir

    store_factory = Callable(
        help="""A callable returning the storage backend for notebook signatures.
         The default uses an SQLite database."""
    ).tag(config=True)

    @default("store_factory")
    def _store_factory_default(self):
        def factory():
            if sqlite3 is None:
                self.log.warning(  # type:ignore[unreachable]
                    "Missing SQLite3, all notebooks will be untrusted!"
                )
                return MemorySignatureStore()
            return SQLiteSignatureStore(self.db_file)

        return factory

    db_file = Unicode(
        help="""The sqlite file in which to store notebook signatures.
        By default, this will be in your Jupyter data directory.
        You can set it to ':memory:' to disable sqlite writing to the filesystem.
        """
    ).tag(config=True)

    @default("db_file")
    def _db_file_default(self):
        if not self.data_dir:
            return ":memory:"
        return str(Path(self.data_dir) / "nbsignatures.db")

    algorithm = Enum(
        algorithms,
        default_value="sha256",
        help="""The hashing algorithm used to sign notebooks.""",
    ).tag(config=True)

    @observe("algorithm")
    def _algorithm_changed(self, change):
        self.digestmod = getattr(hashlib, change["new"])

    digestmod = Any()

    @default("digestmod")
    def _digestmod_default(self):
        return getattr(hashlib, self.algorithm)

    secret_file = Unicode(help="""The file where the secret key is stored.""").tag(config=True)

    @default("secret_file")
    def _secret_file_default(self):
        if not self.data_dir:
            return ""
        return str(Path(self.data_dir) / "notebook_secret")

    secret = Bytes(help="""The secret key with which notebooks are signed.""").tag(config=True)

    @default("secret")
    def _secret_default(self):
        # note : this assumes an Application is running
        if Path(self.secret_file).exists():
            with Path(self.secret_file).open("rb") as f:
                return f.read()
        else:
            secret = encodebytes(os.urandom(1024))
            self._write_secret_file(secret)
            return secret

    def __init__(self, **kwargs):
        """Initialize the notary."""
        super().__init__(**kwargs)
        self.store = self.store_factory()

    def _write_secret_file(self, secret):
        """write my secret to my secret_file"""
        self.log.info("Writing notebook-signing key to %s", self.secret_file)
        with Path(self.secret_file).open("wb") as f:
            f.write(secret)
        try:
            Path(self.secret_file).chmod(0o600)
        except OSError:
            self.log.warning("Could not set permissions on %s", self.secret_file)
        return secret

    def compute_signature(self, nb):
        """Compute a notebook's signature

        by hashing the entire contents of the notebook via HMAC digest.
        """
        hmac = HMAC(self.secret, digestmod=self.digestmod)
        # don't include the previous hash in the content to hash
        with signature_removed(nb):
            # sign the whole thing
            for b in yield_everything(nb):
                hmac.update(b)

        return hmac.hexdigest()

    def check_signature(self, nb):
        """Check a notebook's stored signature

        If a signature is stored in the notebook's metadata,
        a new signature is computed and compared with the stored value.

        Returns True if the signature is found and matches, False otherwise.

        The following conditions must all be met for a notebook to be trusted:
        - a signature is stored in the form 'scheme:hexdigest'
        - the stored scheme matches the requested scheme
        - the requested scheme is available from hashlib
        - the computed hash from notebook_signature matches the stored hash
        """
        if nb.nbformat < 3:
            return False
        signature = self.compute_signature(nb)
        return self.store.check_signature(signature, self.algorithm)

    def sign(self, nb):
        """Sign a notebook, indicating that its output is trusted on this machine

        Stores hash algorithm and hmac digest in a local database of trusted notebooks.
        """
        if nb.nbformat < 3:
            return
        signature = self.compute_signature(nb)
        self.store.store_signature(signature, self.algorithm)

    def unsign(self, nb):
        """Ensure that a notebook is untrusted

        by removing its signature from the trusted database, if present.
        """
        signature = self.compute_signature(nb)
        self.store.remove_signature(signature, self.algorithm)

    def mark_cells(self, nb, trusted):
        """Mark cells as trusted if the notebook's signature can be verified

        Sets ``cell.metadata.trusted = True | False`` on all code cells,
        depending on the *trusted* parameter. This will typically be the return
        value from ``self.check_signature(nb)``.

        This function is the inverse of check_cells
        """
        if nb.nbformat < 3:
            return

        for cell in yield_code_cells(nb):
            cell["metadata"]["trusted"] = trusted

    def _check_cell(self, cell, nbformat_version):
        """Do we trust an individual cell?

        Return True if:

        - cell is explicitly trusted
        - cell has no potentially unsafe rich output

        If a cell has no output, or only simple print statements,
        it will always be trusted.
        """
        # explicitly trusted
        if cell["metadata"].pop("trusted", False):
            return True

        # explicitly safe output
        if nbformat_version >= 4:
            unsafe_output_types = ["execute_result", "display_data"]
            safe_keys = {"output_type", "execution_count", "metadata"}
        else:  # v3
            unsafe_output_types = ["pyout", "display_data"]
            safe_keys = {"output_type", "prompt_number", "metadata"}

        for output in cell["outputs"]:
            output_type = output["output_type"]
            if output_type in unsafe_output_types:
                # if there are any data keys not in the safe whitelist
                output_keys = set(output)
                if output_keys.difference(safe_keys):
                    return False

        return True

    def check_cells(self, nb):
        """Return whether all code cells are trusted.

        A cell is trusted if the 'trusted' field in its metadata is truthy, or
        if it has no potentially unsafe outputs.
        If there are no code cells, return True.

        This function is the inverse of mark_cells.
        """
        if nb.nbformat < 3:
            return False
        trusted = True
        for cell in yield_code_cells(nb):
            # only distrust a cell if it actually has some output to distrust
            if not self._check_cell(cell, nb.nbformat):
                trusted = False

        return trusted


trust_flags: dict[str, t.Any] = {
    "reset": (
        {"TrustNotebookApp": {"reset": True}},
        """Delete the trusted notebook cache.
        All previously signed notebooks will become untrusted.
        """,
    ),
}
trust_flags.update(base_flags)


class TrustNotebookApp(JupyterApp):
    """An application for handling notebook trust."""

    version = __version__
    description = """Sign one or more Jupyter notebooks with your key,
    to trust their dynamic (HTML, Javascript) output.

    Otherwise, you will have to re-execute the notebook to see output.
    """
    # This command line tool should use the same config file as the notebook

    @default("config_file_name")
    def _config_file_name_default(self):
        return "jupyter_notebook_config"

    examples = """
    jupyter trust mynotebook.ipynb and_this_one.ipynb
    """

    flags = trust_flags

    reset = Bool(
        False,
        help="""If True, delete the trusted signature cache.
        After reset, all previously signed notebooks will become untrusted.
        """,
    ).tag(config=True)

    notary = Instance(NotebookNotary)

    @default("notary")
    def _notary_default(self):
        return NotebookNotary(parent=self, data_dir=self.data_dir)

    def sign_notebook_file(self, notebook_path):
        """Sign a notebook from the filesystem"""
        if not Path(notebook_path).exists():
            self.log.error("Notebook missing: %s", notebook_path)
            self.exit(1)
        with Path(notebook_path).open(encoding="utf8") as f:
            nb = read(f, NO_CONVERT)
        self.sign_notebook(nb, notebook_path)

    def sign_notebook(self, nb, notebook_path="<stdin>"):
        """Sign a notebook that's been loaded"""
        if self.notary.check_signature(nb):
            print("Notebook already signed: %s" % notebook_path)  # noqa: T201
        else:
            print("Signing notebook: %s" % notebook_path)  # noqa: T201
            self.notary.sign(nb)

    def generate_new_key(self):
        """Generate a new notebook signature key"""
        print("Generating new notebook key: %s" % self.notary.secret_file)  # noqa: T201
        self.notary._write_secret_file(os.urandom(1024))

    def start(self):
        """Start the trust notebook app."""
        if self.reset:
            if Path(self.notary.db_file).exists():
                print("Removing trusted signature cache: %s" % self.notary.db_file)  # noqa: T201
                Path(self.notary.db_file).unlink()
            self.generate_new_key()
            return
        if not self.extra_args:
            self.log.debug("Reading notebook from stdin")
            nb_s = sys.stdin.read()
            assert isinstance(nb_s, str)
            nb = reads(nb_s, NO_CONVERT)
            self.sign_notebook(nb, "<stdin>")
        else:
            for notebook_path in self.extra_args:
                self.sign_notebook_file(notebook_path)


main = TrustNotebookApp.launch_instance

if __name__ == "__main__":
    main()
