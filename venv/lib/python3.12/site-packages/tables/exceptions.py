"""Declare exceptions and warnings that are specific to PyTables."""

from __future__ import annotations

import os
import warnings
import traceback
from collections.abc import Callable

__all__ = [
    "ChunkError",
    "ClosedFileError",
    "ClosedNodeError",
    "DataTypeWarning",
    "ExperimentalFeatureWarning",
    "FileModeError",
    "FiltersWarning",
    "FlavorError",
    "FlavorWarning",
    "HDF5ExtError",
    "NaturalNameWarning",
    "NoSuchChunkError",
    "NoSuchNodeError",
    "NodeError",
    "NotChunkedError",
    "NotChunkAlignedError",
    "OldIndexWarning",
    "PerformanceWarning",
    "UnclosedFileWarning",
    "UndoRedoError",
    "UndoRedoWarning",
]


__docformat__ = "reStructuredText"
"""The format of documentation strings in this module."""


class HDF5ExtError(RuntimeError):
    """A low level HDF5 operation failed.

    This exception is raised the low level PyTables components used for
    accessing HDF5 files.  It usually signals that something is not
    going well in the HDF5 library or even at the Input/Output level.

    Errors in the HDF5 C library may be accompanied by an extensive
    HDF5 back trace on standard error (see also
    :func:`tables.silence_hdf5_messages`).

    .. versionchanged:: 2.4

    Parameters
    ----------
    message
        error message
    h5bt
        This parameter (keyword only) controls the HDF5 back trace
        handling. Any keyword arguments other than h5bt is ignored.

        * if set to False the HDF5 back trace is ignored and the
          :attr:`HDF5ExtError.h5backtrace` attribute is set to None
        * if set to True the back trace is retrieved from the HDF5
          library and stored in the :attr:`HDF5ExtError.h5backtrace`
          attribute as a list of tuples
        * if set to "VERBOSE" (default) the HDF5 back trace is
          stored in the :attr:`HDF5ExtError.h5backtrace` attribute
          and also included in the string representation of the
          exception
        * if not set (or set to None) the default policy is used
          (see :attr:`HDF5ExtError.DEFAULT_H5_BACKTRACE_POLICY`)

    """

    # NOTE: in order to avoid circular dependencies between modules the
    #       _dump_h5_backtrace method is set at initialization time in
    #       the utilsextension.pyx.
    _dump_h5_backtrace: (
        Callable[[], list[tuple[str, int, str, str]]] | None
    ) = None

    DEFAULT_H5_BACKTRACE_POLICY = "VERBOSE"
    """Default policy for HDF5 backtrace handling

    * if set to False the HDF5 back trace is ignored and the
      :attr:`HDF5ExtError.h5backtrace` attribute is set to None
    * if set to True the back trace is retrieved from the HDF5
      library and stored in the :attr:`HDF5ExtError.h5backtrace`
      attribute as a list of tuples
    * if set to "VERBOSE" (default) the HDF5 back trace is
      stored in the :attr:`HDF5ExtError.h5backtrace` attribute
      and also included in the string representation of the
      exception

    This parameter can be set using the
    :envvar:`PT_DEFAULT_H5_BACKTRACE_POLICY` environment variable.
    Allowed values are "IGNORE" (or "FALSE"), "SAVE" (or "TRUE") and
    "VERBOSE" to set the policy to False, True and "VERBOSE"
    respectively.  The special value "DEFAULT" can be used to reset
    the policy to the default value

    .. versionadded:: 2.4
    """

    @classmethod
    def set_policy_from_env(cls) -> str:
        """Set the policy from environment variables."""
        envmap = {
            "IGNORE": False,
            "FALSE": False,
            "SAVE": True,
            "TRUE": True,
            "VERBOSE": "VERBOSE",
            "DEFAULT": "VERBOSE",
        }
        oldvalue = cls.DEFAULT_H5_BACKTRACE_POLICY
        envvalue = os.environ.get("PT_DEFAULT_H5_BACKTRACE_POLICY", "DEFAULT")
        try:
            newvalue = envmap[envvalue.upper()]
        except KeyError:
            warnings.warn(
                "Invalid value for the environment variable "
                "'PT_DEFAULT_H5_BACKTRACE_POLICY'.  The default "
                "policy for HDF5 back trace management in PyTables "
                "will be: '%s'" % oldvalue
            )
        else:
            cls.DEFAULT_H5_BACKTRACE_POLICY = newvalue

        return oldvalue

    def __init__(self, *args, **kargs) -> None:

        super().__init__(*args)

        self._h5bt_policy = kargs.get("h5bt", self.DEFAULT_H5_BACKTRACE_POLICY)

        if self._h5bt_policy and self._dump_h5_backtrace is not None:
            self.h5backtrace = self._dump_h5_backtrace()
            """HDF5 back trace.

            Contains the HDF5 back trace as a (possibly empty) list of
            tuples.  Each tuple has the following format::

                (filename, line number, function name, text)

            Depending on the value of the *h5bt* parameter passed to the
            initializer the h5backtrace attribute can be set to None.
            This means that the HDF5 back trace has been simply ignored
            (not retrieved from the HDF5 C library error stack) or that
            there has been an error (silently ignored) during the HDF5 back
            trace retrieval.

            .. versionadded:: 2.4

            See Also
            --------
            traceback.format_list : :func:`traceback.format_list`

            """

            # XXX: check _dump_h5_backtrace failures
        else:
            self.h5backtrace = None

    def __str__(self) -> str:
        """Return a sting representation of the exception.

        The actual result depends on policy set in the initializer
        :meth:`HDF5ExtError.__init__`.

        .. versionadded:: 2.4

        """
        verbose = bool(self._h5bt_policy in ("VERBOSE", "verbose"))

        if verbose and self.h5backtrace:
            bt = "\n".join(
                [
                    "HDF5 error back trace\n",
                    self.format_h5_backtrace(),
                    "End of HDF5 error back trace",
                ]
            )

            if len(self.args) == 1 and isinstance(self.args[0], str):
                msg = super().__str__()
                msg = f"{bt}\n\n{msg}"
            elif self.h5backtrace[-1][-1]:
                msg = f"{bt}\n\n{self.h5backtrace[-1][-1]}"
            else:
                msg = bt
        else:
            msg = super().__str__()

        return msg

    def format_h5_backtrace(
        self, backtrace: list[tuple[str, int, str, str]] | None = None
    ) -> str:
        """Convert the HDF5 trace back into a string.

        The HDF5 trace back is represented as a list of tuples.

        See :attr:`HDF5ExtError.h5backtrace`.

        .. versionadded:: 2.4

        """
        if backtrace is None:
            backtrace = self.h5backtrace

        if backtrace is None:
            return "No HDF5 back trace available"
        else:
            return "".join(traceback.format_list(backtrace))


# Initialize the policy for HDF5 back trace handling
HDF5ExtError.set_policy_from_env()


# The following exceptions are concretions of the ``ValueError`` exceptions
# raised by ``file`` objects on certain operations.


class ClosedNodeError(ValueError):
    """The operation can not be completed because the node is closed.

    For instance, listing the children of a closed group is not allowed.

    """

    pass


class ClosedFileError(ValueError):
    """The operation can not be completed because the hosting file is closed.

    For instance, getting an existing node from a closed file is not
    allowed.

    """

    pass


class FileModeError(ValueError):
    """FIle mode error.

    The operation can not be carried out because the mode in which the
    hosting file is opened is not adequate.

    For instance, removing an existing leaf from a read-only file is not
    allowed.

    """

    pass


class NodeError(AttributeError, LookupError):
    """Invalid hierarchy manipulation operation requested.

    This exception is raised when the user requests an operation on the
    hierarchy which can not be run because of the current layout of the
    tree.  This includes accessing nonexistent nodes, moving or copying
    or creating over an existing node, non-recursively removing groups
    with children, and other similarly invalid operations.

    A node in a PyTables database cannot be simply overwritten by
    replacing it.  Instead, the old node must be removed explicitly
    before another one can take its place.  This is done to protect
    interactive users from inadvertently deleting whole trees of data by
    a single erroneous command.

    """

    pass


class NoSuchNodeError(NodeError):
    """An operation was requested on a node that does not exist.

    This exception is raised when an operation gets a path name or a
    ``(where, name)`` pair leading to a nonexistent node.

    """

    pass


class UndoRedoError(Exception):
    """Problems with doing/redoing actions with Undo/Redo feature.

    This exception indicates a problem related to the Undo/Redo
    mechanism, such as trying to undo or redo actions with this
    mechanism disabled, or going to a nonexistent mark.

    """

    pass


class UndoRedoWarning(Warning):
    """Issued when an action not supporting Undo/Redo is run.

    This warning is only shown when the Undo/Redo mechanism is enabled.

    """

    pass


class NaturalNameWarning(Warning):
    """Issued when a non-pythonic name is given for a node.

    This is not an error and may even be very useful in certain
    contexts, but one should be aware that such nodes cannot be
    accessed using natural naming (instead, ``getattr()`` must be
    used explicitly).
    """

    pass


class PerformanceWarning(Warning):
    """Warning for operations which may cause a performance drop.

    This warning is issued when an operation is made on the database
    which may cause it to slow down on future operations (i.e. making
    the node tree grow too much).

    """

    pass


class FlavorError(ValueError):
    """Unsupported or unavailable flavor or flavor conversion.

    This exception is raised when an unsupported or unavailable flavor
    is given to a dataset, or when a conversion of data between two
    given flavors is not supported nor available.

    """

    pass


class FlavorWarning(Warning):
    """Unsupported or unavailable flavor conversion.

    This warning is issued when a conversion of data between two given
    flavors is not supported nor available, and raising an error would
    render the data inaccessible (e.g. on a dataset of an unavailable
    flavor in a read-only file).

    See the `FlavorError` class for more information.

    """

    pass


class FiltersWarning(Warning):
    """Unavailable filters.

    This warning is issued when a valid filter is specified but it is
    not available in the system.  It may mean that an available default
    filter is to be used instead.

    """

    pass


class OldIndexWarning(Warning):
    """Unsupported index format.

    This warning is issued when an index in an unsupported format is
    found.  The index will be marked as invalid and will behave as if
    it doesn't exist.

    """

    pass


class DataTypeWarning(Warning):
    """Unsupported data type.

    This warning is issued when an unsupported HDF5 data type is found
    (normally in a file created with other tool than PyTables).

    """

    pass


class ExperimentalFeatureWarning(Warning):
    """Generic warning for experimental features.

    This warning is issued when using a functionality that is still
    experimental and that users have to use with care.

    """

    pass


class UnclosedFileWarning(Warning):
    """Warning raised when there are still open files at program exit.

    PyTables will close remaining open files at exit, but raise this warning.
    """

    pass


class ChunkError(Exception):
    """An operation related to direct chunk access failed.

    This exception may be related with the properties of the dataset or the
    chunk being accessed, or with how the chunk is being accessed.  It is a
    base for more specific exceptions.

    """

    pass


class NotChunkedError(ChunkError):
    """A direct chunking operation was attempted on a non-chunked dataset.

    For instance, chunk information was requested for a plain ``Array``
    instance.

    """

    pass


class NotChunkAlignedError(ChunkError):
    """Coordinate not aligned to the chunks.

    A direct chunk read/write operation was given coordinates that do not
    match the chunk's start.

    These operations require coordinates that are integer multiples of the
    dataset's chunksize.

    """

    pass


class NoSuchChunkError(ChunkError):
    """The chunk with the given coordinates does not exist in storage.

    The coordinates are within the dataset's shape, though.

    This is only an error when the chunk is to be read.  Such a missing chunk
    can be written, in which case it is created in storage.

    """

    pass
