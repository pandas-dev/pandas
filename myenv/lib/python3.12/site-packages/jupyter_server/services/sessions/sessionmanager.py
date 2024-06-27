"""A base class session manager."""

# Copyright (c) Jupyter Development Team.
# Distributed under the terms of the Modified BSD License.
import os
import pathlib
import uuid
from typing import Any, Dict, List, NewType, Optional, Union, cast

KernelName = NewType("KernelName", str)
ModelName = NewType("ModelName", str)

try:
    import sqlite3
except ImportError:
    # fallback on pysqlite2 if Python was build without sqlite
    from pysqlite2 import dbapi2 as sqlite3  # type:ignore[no-redef]

from dataclasses import dataclass, fields

from jupyter_core.utils import ensure_async
from tornado import web
from traitlets import Instance, TraitError, Unicode, validate
from traitlets.config.configurable import LoggingConfigurable

from jupyter_server.traittypes import InstanceFromClasses


class KernelSessionRecordConflict(Exception):
    """Exception class to use when two KernelSessionRecords cannot
    merge because of conflicting data.
    """


@dataclass
class KernelSessionRecord:
    """A record object for tracking a Jupyter Server Kernel Session.

    Two records that share a session_id must also share a kernel_id, while
    kernels can have multiple session (and thereby) session_ids
    associated with them.
    """

    session_id: Optional[str] = None
    kernel_id: Optional[str] = None

    def __eq__(self, other: object) -> bool:
        """Whether a record equals another."""
        if isinstance(other, KernelSessionRecord):
            condition1 = self.kernel_id and self.kernel_id == other.kernel_id
            condition2 = all(
                [
                    self.session_id == other.session_id,
                    self.kernel_id is None or other.kernel_id is None,
                ]
            )
            if any([condition1, condition2]):
                return True
            # If two records share session_id but have different kernels, this is
            # and ill-posed expression. This should never be true. Raise an exception
            # to inform the user.
            if all(
                [
                    self.session_id,
                    self.session_id == other.session_id,
                    self.kernel_id != other.kernel_id,
                ]
            ):
                msg = (
                    "A single session_id can only have one kernel_id "
                    "associated with. These two KernelSessionRecords share the same "
                    "session_id but have different kernel_ids. This should "
                    "not be possible and is likely an issue with the session "
                    "records."
                )
                raise KernelSessionRecordConflict(msg)
        return False

    def update(self, other: "KernelSessionRecord") -> None:
        """Updates in-place a kernel from other (only accepts positive updates"""
        if not isinstance(other, KernelSessionRecord):
            msg = "'other' must be an instance of KernelSessionRecord."  # type:ignore[unreachable]
            raise TypeError(msg)

        if other.kernel_id and self.kernel_id and other.kernel_id != self.kernel_id:
            msg = "Could not update the record from 'other' because the two records conflict."
            raise KernelSessionRecordConflict(msg)

        for field in fields(self):
            if hasattr(other, field.name) and getattr(other, field.name):
                setattr(self, field.name, getattr(other, field.name))


class KernelSessionRecordList:
    """An object for storing and managing a list of KernelSessionRecords.

    When adding a record to the list, the KernelSessionRecordList
    first checks if the record already exists in the list. If it does,
    the record will be updated with the new information; otherwise,
    it will be appended.
    """

    _records: List[KernelSessionRecord]

    def __init__(self, *records: KernelSessionRecord):
        """Initialize a record list."""
        self._records = []
        for record in records:
            self.update(record)

    def __str__(self):
        """The string representation of a record list."""
        return str(self._records)

    def __contains__(self, record: Union[KernelSessionRecord, str]) -> bool:
        """Search for records by kernel_id and session_id"""
        if isinstance(record, KernelSessionRecord) and record in self._records:
            return True

        if isinstance(record, str):
            for r in self._records:
                if record in [r.session_id, r.kernel_id]:
                    return True
        return False

    def __len__(self):
        """The length of the record list."""
        return len(self._records)

    def get(self, record: Union[KernelSessionRecord, str]) -> KernelSessionRecord:
        """Return a full KernelSessionRecord from a session_id, kernel_id, or
        incomplete KernelSessionRecord.
        """
        if isinstance(record, str):
            for r in self._records:
                if record in (r.kernel_id, r.session_id):
                    return r
        elif isinstance(record, KernelSessionRecord):
            for r in self._records:
                if record == r:
                    return record
        msg = f"{record} not found in KernelSessionRecordList."
        raise ValueError(msg)

    def update(self, record: KernelSessionRecord) -> None:
        """Update a record in-place or append it if not in the list."""
        try:
            idx = self._records.index(record)
            self._records[idx].update(record)
        except ValueError:
            self._records.append(record)

    def remove(self, record: KernelSessionRecord) -> None:
        """Remove a record if its found in the list. If it's not found,
        do nothing.
        """
        if record in self._records:
            self._records.remove(record)


class SessionManager(LoggingConfigurable):
    """A session manager."""

    database_filepath = Unicode(
        default_value=":memory:",
        help=(
            "The filesystem path to SQLite Database file "
            "(e.g. /path/to/session_database.db). By default, the session "
            "database is stored in-memory (i.e. `:memory:` setting from sqlite3) "
            "and does not persist when the current Jupyter Server shuts down."
        ),
    ).tag(config=True)

    @validate("database_filepath")
    def _validate_database_filepath(self, proposal):
        """Validate a database file path."""
        value = proposal["value"]
        if value == ":memory:":
            return value
        path = pathlib.Path(value)
        if path.exists():
            # Verify that the database path is not a directory.
            if path.is_dir():
                msg = "`database_filepath` expected a file path, but the given path is a directory."
                raise TraitError(msg)
            # Verify that database path is an SQLite 3 Database by checking its header.
            with open(value, "rb") as f:
                header = f.read(100)

            if not header.startswith(b"SQLite format 3") and header != b"":
                msg = "The given file is not an SQLite database file."
                raise TraitError(msg)
        return value

    kernel_manager = Instance("jupyter_server.services.kernels.kernelmanager.MappingKernelManager")
    contents_manager = InstanceFromClasses(
        [
            "jupyter_server.services.contents.manager.ContentsManager",
            "notebook.services.contents.manager.ContentsManager",
        ]
    )

    def __init__(self, *args, **kwargs):
        """Initialize a record list."""
        super().__init__(*args, **kwargs)
        self._pending_sessions = KernelSessionRecordList()

    # Session database initialized below
    _cursor = None
    _connection = None
    _columns = {"session_id", "path", "name", "type", "kernel_id"}

    @property
    def cursor(self):
        """Start a cursor and create a database called 'session'"""
        if self._cursor is None:
            self._cursor = self.connection.cursor()
            self._cursor.execute(
                """CREATE TABLE IF NOT EXISTS session
                (session_id, path, name, type, kernel_id)"""
            )
        return self._cursor

    @property
    def connection(self):
        """Start a database connection"""
        if self._connection is None:
            # Set isolation level to None to autocommit all changes to the database.
            self._connection = sqlite3.connect(self.database_filepath, isolation_level=None)
            self._connection.row_factory = sqlite3.Row
        return self._connection

    def close(self):
        """Close the sqlite connection"""
        if self._cursor is not None:
            self._cursor.close()
            self._cursor = None

    def __del__(self):
        """Close connection once SessionManager closes"""
        self.close()

    async def session_exists(self, path):
        """Check to see if the session of a given name exists"""
        exists = False
        self.cursor.execute("SELECT * FROM session WHERE path=?", (path,))
        row = self.cursor.fetchone()
        if row is not None:
            # Note, although we found a row for the session, the associated kernel may have
            # been culled or died unexpectedly.  If that's the case, we should delete the
            # row, thereby terminating the session.  This can be done via a call to
            # row_to_model that tolerates that condition.  If row_to_model returns None,
            # we'll return false, since, at that point, the session doesn't exist anyway.
            model = await self.row_to_model(row, tolerate_culled=True)
            if model is not None:
                exists = True
        return exists

    def new_session_id(self) -> str:
        """Create a uuid for a new session"""
        return str(uuid.uuid4())

    async def create_session(
        self,
        path: Optional[str] = None,
        name: Optional[ModelName] = None,
        type: Optional[str] = None,
        kernel_name: Optional[KernelName] = None,
        kernel_id: Optional[str] = None,
    ) -> Dict[str, Any]:
        """Creates a session and returns its model

        Parameters
        ----------
        name: ModelName(str)
            Usually the model name, like the filename associated with current
            kernel.
        """
        session_id = self.new_session_id()
        record = KernelSessionRecord(session_id=session_id)
        self._pending_sessions.update(record)
        if kernel_id is not None and kernel_id in self.kernel_manager:
            pass
        else:
            kernel_id = await self.start_kernel_for_session(
                session_id, path, name, type, kernel_name
            )
        record.kernel_id = kernel_id
        self._pending_sessions.update(record)
        result = await self.save_session(
            session_id, path=path, name=name, type=type, kernel_id=kernel_id
        )
        self._pending_sessions.remove(record)
        return cast(Dict[str, Any], result)

    def get_kernel_env(
        self, path: Optional[str], name: Optional[ModelName] = None
    ) -> Dict[str, str]:
        """Return the environment variables that need to be set in the kernel

        Parameters
        ----------
        path : str
            the url path for the given session.
        name: ModelName(str), optional
            Here the name is likely to be the name of the associated file
            with the current kernel at startup time.
        """
        if name is not None:
            cwd = self.kernel_manager.cwd_for_path(path)
            path = os.path.join(cwd, name)
        assert isinstance(path, str)
        return {**os.environ, "JPY_SESSION_NAME": path}

    async def start_kernel_for_session(
        self,
        session_id: str,
        path: Optional[str],
        name: Optional[ModelName],
        type: Optional[str],
        kernel_name: Optional[KernelName],
    ) -> str:
        """Start a new kernel for a given session.

        Parameters
        ----------
        session_id : str
            uuid for the session; this method must be given a session_id
        path : str
            the path for the given session - seem to be a session id sometime.
        name : str
            Usually the model name, like the filename associated with current
            kernel.
        type : str
            the type of the session
        kernel_name : str
            the name of the kernel specification to use.  The default kernel name will be used if not provided.
        """
        # allow contents manager to specify kernels cwd
        kernel_path = await ensure_async(self.contents_manager.get_kernel_path(path=path))

        kernel_env = self.get_kernel_env(path, name)
        kernel_id = await self.kernel_manager.start_kernel(
            path=kernel_path,
            kernel_name=kernel_name,
            env=kernel_env,
        )
        return cast(str, kernel_id)

    async def save_session(self, session_id, path=None, name=None, type=None, kernel_id=None):
        """Saves the items for the session with the given session_id

        Given a session_id (and any other of the arguments), this method
        creates a row in the sqlite session database that holds the information
        for a session.

        Parameters
        ----------
        session_id : str
            uuid for the session; this method must be given a session_id
        path : str
            the path for the given session
        name : str
            the name of the session
        type : str
            the type of the session
        kernel_id : str
            a uuid for the kernel associated with this session

        Returns
        -------
        model : dict
            a dictionary of the session model
        """
        self.cursor.execute(
            "INSERT INTO session VALUES (?,?,?,?,?)",
            (session_id, path, name, type, kernel_id),
        )
        result = await self.get_session(session_id=session_id)
        return result

    async def get_session(self, **kwargs):
        """Returns the model for a particular session.

        Takes a keyword argument and searches for the value in the session
        database, then returns the rest of the session's info.

        Parameters
        ----------
        **kwargs : dict
            must be given one of the keywords and values from the session database
            (i.e. session_id, path, name, type, kernel_id)

        Returns
        -------
        model : dict
            returns a dictionary that includes all the information from the
            session described by the kwarg.
        """
        if not kwargs:
            msg = "must specify a column to query"
            raise TypeError(msg)

        conditions = []
        for column in kwargs:
            if column not in self._columns:
                msg = f"No such column: {column}"
                raise TypeError(msg)
            conditions.append("%s=?" % column)

        query = "SELECT * FROM session WHERE %s" % (" AND ".join(conditions))  # noqa: S608

        self.cursor.execute(query, list(kwargs.values()))
        try:
            row = self.cursor.fetchone()
        except KeyError:
            # The kernel is missing, so the session just got deleted.
            row = None

        if row is None:
            q = []
            for key, value in kwargs.items():
                q.append(f"{key}={value!r}")

            raise web.HTTPError(404, "Session not found: %s" % (", ".join(q)))

        try:
            model = await self.row_to_model(row)
        except KeyError as e:
            raise web.HTTPError(404, "Session not found: %s" % str(e)) from e
        return model

    async def update_session(self, session_id, **kwargs):
        """Updates the values in the session database.

        Changes the values of the session with the given session_id
        with the values from the keyword arguments.

        Parameters
        ----------
        session_id : str
            a uuid that identifies a session in the sqlite3 database
        **kwargs : str
            the key must correspond to a column title in session database,
            and the value replaces the current value in the session
            with session_id.
        """
        await self.get_session(session_id=session_id)

        if not kwargs:
            # no changes
            return

        sets = []
        for column in kwargs:
            if column not in self._columns:
                raise TypeError("No such column: %r" % column)
            sets.append("%s=?" % column)
        query = "UPDATE session SET %s WHERE session_id=?" % (", ".join(sets))  # noqa: S608
        self.cursor.execute(query, [*list(kwargs.values()), session_id])

        if hasattr(self.kernel_manager, "update_env"):
            self.cursor.execute(
                "SELECT path, name, kernel_id FROM session WHERE session_id=?", [session_id]
            )
            path, name, kernel_id = self.cursor.fetchone()
            self.kernel_manager.update_env(kernel_id=kernel_id, env=self.get_kernel_env(path, name))

    async def kernel_culled(self, kernel_id: str) -> bool:
        """Checks if the kernel is still considered alive and returns true if its not found."""
        return kernel_id not in self.kernel_manager

    async def row_to_model(self, row, tolerate_culled=False):
        """Takes sqlite database session row and turns it into a dictionary"""
        kernel_culled: bool = await ensure_async(self.kernel_culled(row["kernel_id"]))
        if kernel_culled:
            # The kernel was culled or died without deleting the session.
            # We can't use delete_session here because that tries to find
            # and shut down the kernel - so we'll delete the row directly.
            #
            # If caller wishes to tolerate culled kernels, log a warning
            # and return None.  Otherwise, raise KeyError with a similar
            # message.
            self.cursor.execute("DELETE FROM session WHERE session_id=?", (row["session_id"],))
            msg = (
                "Kernel '{kernel_id}' appears to have been culled or died unexpectedly, "
                "invalidating session '{session_id}'. The session has been removed.".format(
                    kernel_id=row["kernel_id"], session_id=row["session_id"]
                )
            )
            if tolerate_culled:
                self.log.warning(f"{msg}  Continuing...")
                return None
            raise KeyError(msg)

        kernel_model = await ensure_async(self.kernel_manager.kernel_model(row["kernel_id"]))
        model = {
            "id": row["session_id"],
            "path": row["path"],
            "name": row["name"],
            "type": row["type"],
            "kernel": kernel_model,
        }
        if row["type"] == "notebook":
            # Provide the deprecated API.
            model["notebook"] = {"path": row["path"], "name": row["name"]}
        return model

    async def list_sessions(self):
        """Returns a list of dictionaries containing all the information from
        the session database"""
        c = self.cursor.execute("SELECT * FROM session")
        result = []
        # We need to use fetchall() here, because row_to_model can delete rows,
        # which messes up the cursor if we're iterating over rows.
        for row in c.fetchall():
            try:
                model = await self.row_to_model(row)
                result.append(model)
            except KeyError:
                pass
        return result

    async def delete_session(self, session_id):
        """Deletes the row in the session database with given session_id"""
        record = KernelSessionRecord(session_id=session_id)
        self._pending_sessions.update(record)
        session = await self.get_session(session_id=session_id)
        await ensure_async(self.kernel_manager.shutdown_kernel(session["kernel"]["id"]))
        self.cursor.execute("DELETE FROM session WHERE session_id=?", (session_id,))
        self._pending_sessions.remove(record)
