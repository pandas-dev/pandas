"""Parameters for PyTables."""

import os as _os

__docformat__ = "reStructuredText"
"""The format of documentation strings in this module."""

_KB = 1024
"""The size of a Kilobyte in bytes"""

_MB = 1024 * _KB
"""The size of a Megabyte in bytes"""

# Tunable parameters
# ==================
# Be careful when touching these!

# Parameters for different internal caches
# ----------------------------------------

BOUNDS_MAX_SIZE = 1 * _MB
"""The maximum size for bounds values cached during index lookups."""

BOUNDS_MAX_SLOTS = 4 * _KB
"""The maximum number of slots for the BOUNDS cache."""

ITERSEQ_MAX_ELEMENTS = 1 * _KB
"""The maximum number of iterator elements cached in data lookups."""

ITERSEQ_MAX_SIZE = 1 * _MB
"""The maximum space that will take ITERSEQ cache (in bytes)."""

ITERSEQ_MAX_SLOTS = 128
"""The maximum number of slots in ITERSEQ cache."""

LIMBOUNDS_MAX_SIZE = 256 * _KB
"""The maximum size for the query limits (for example, ``(lim1, lim2)``
in conditions like ``lim1 <= col < lim2``) cached during index lookups
(in bytes)."""

LIMBOUNDS_MAX_SLOTS = 128
"""The maximum number of slots for LIMBOUNDS cache."""

TABLE_MAX_SIZE = 1 * _MB
"""The maximum size for table chunks cached during index queries."""

SORTED_MAX_SIZE = 1 * _MB
"""The maximum size for sorted values cached during index lookups."""

SORTEDLR_MAX_SIZE = 8 * _MB
"""The maximum size for chunks in last row cached in index lookups (in
bytes)."""

SORTEDLR_MAX_SLOTS = 1 * _KB
"""The maximum number of chunks for SORTEDLR cache."""


# Parameters for general cache behaviour
# --------------------------------------
#
# The next parameters will not be effective if passed to the
# `open_file()` function (so, they can only be changed in a *global*
# way).  You can change them in the file, but this is strongly
# discouraged unless you know well what you are doing.

DISABLE_EVERY_CYCLES = 10
"""The number of cycles in which a cache will be forced to be disabled
if the hit ratio is lower than the LOWEST_HIT_RATIO (see below).  This
value should provide time enough to check whether the cache is being
efficient or not."""

ENABLE_EVERY_CYCLES = 50
"""The number of cycles in which a cache will be forced to be
(re-)enabled, irregardless of the hit ratio. This will provide a chance
for checking if we are in a better scenario for doing caching again."""

LOWEST_HIT_RATIO = 0.6
"""The minimum acceptable hit ratio for a cache to avoid disabling (and
freeing) it."""


# Tunable parameters
# ==================
# Be careful when touching these!

# Recommended maximum values
# --------------------------

# Following are the recommended values for several limits.  However,
# these limits are somewhat arbitrary and can be increased if you have
# enough resources.

MAX_COLUMNS = 512
"""Maximum number of columns in :class:`tables.Table` objects before a
:exc:`tables.PerformanceWarning` is issued.  This limit is somewhat
arbitrary and can be increased.
"""

MAX_NODE_ATTRS = 4 * _KB
"""Maximum allowed number of attributes in a node."""

MAX_GROUP_WIDTH = 16 * _KB
"""Maximum allowed number of children hanging from a group."""

MAX_TREE_DEPTH = 2 * _KB
"""Maximum depth in object tree allowed."""

MAX_UNDO_PATH_LENGTH = 10 * _KB
"""Maximum length of paths allowed in undo/redo operations."""


# Cache limits
# ------------

COND_CACHE_SLOTS = 128
"""Maximum number of conditions for table queries to be kept in memory."""

CHUNK_CACHE_NELMTS = 521
"""Number of elements for HDF5 chunk cache."""

CHUNK_CACHE_PREEMPT = 0.0
"""Chunk preemption policy.  This value should be between 0 and 1
inclusive and indicates how much chunks that have been fully read are
favored for preemption. A value of zero means fully read chunks are
treated no differently than other chunks (the preemption is strictly
LRU) while a value of one means fully read chunks are always preempted
before other chunks."""

CHUNK_CACHE_SIZE = 16 * _MB
"""Size (in bytes) for HDF5 chunk cache."""

# Size for new metadata cache system
METADATA_CACHE_SIZE = 1 * _MB  # 1 MB is the default for HDF5
"""Size (in bytes) of the HDF5 metadata cache."""


# NODE_CACHE_SLOTS tells the number of nodes that fits in the cache.
#
# There are several forces driving the election of this number:
# 1.- As more nodes, better chances to re-use nodes
#     --> better performance
# 2.- As more nodes, the re-ordering of the LRU cache takes more time
#     --> less performance
# 3.- As more nodes, the memory needs for PyTables grows, specially for table
#     writings (that could take double of memory than table reads!).
#
# The default value here is quite conservative. If you have a system
# with tons of memory, and if you are touching regularly a very large
# number of leaves, try increasing this value and see if it fits better
# for you. Please report back your feedback.
# Using a lower value than 64 provides a workaround for issue #977.
NODE_CACHE_SLOTS = 32
"""Maximum number of nodes to be kept in the metadata cache.

It is the number of nodes to be kept in the metadata cache. Least recently
used nodes are unloaded from memory when this number of loaded nodes is
reached. To load a node again, simply access it as usual.
Nodes referenced by user variables and, in general, all nodes that are still
open are registered in the node manager and can be quickly accessed even
if they are not in the cache.

Negative value means that all the touched nodes will be kept in an
internal dictionary.  This is the faster way to load/retrieve nodes.
However, and in order to avoid a large memory consumption, the user will
be warned when the number of loaded nodes will reach the
``-NODE_CACHE_SLOTS`` value.

Finally, a value of zero means that any cache mechanism is disabled.
"""


# Parameters for the I/O buffer in `Leaf` objects
# -----------------------------------------------

IO_BUFFER_SIZE = 16 * _MB
"""The PyTables internal buffer size for I/O purposes.  Should not
exceed the amount of highest level cache size in your CPU."""

BUFFER_TIMES = 100
"""The maximum buffersize/rowsize ratio before issuing a
:exc:`tables.PerformanceWarning`."""


# Miscellaneous
# -------------

EXPECTED_ROWS_EARRAY = 1000
"""Default expected number of rows for :class:`EArray` objects."""

EXPECTED_ROWS_VLARRAY = 1000
"""Default expected number of rows for :class:`VLArray` objects.

.. versionadded:: 3.0

"""

EXPECTED_ROWS_TABLE = 10_000
"""Default expected number of rows for :class:`Table` objects."""

PYTABLES_SYS_ATTRS = True
"""Set this to ``False`` if you don't want to create PyTables system
attributes in datasets.  Also, if set to ``False`` the possible existing
system attributes are not considered for guessing the class of the node
during its loading from disk (this work is delegated to the PyTables'
class discoverer function for general HDF5 files)."""

MAX_NUMEXPR_THREADS = int(_os.environ.get("NUMEXPR_MAX_THREADS", 4))
"""The maximum number of threads that PyTables should use internally in
Numexpr.  If `None`, it is automatically set to the number of cores in
your machine. In general, it is a good idea to set this to the number of
cores in your machine or, when your machine has many of them (e.g. > 8),
perhaps stay at 8 at maximum.  In general, 4 threads is a good tradeoff."""

MAX_BLOSC_THREADS = 1  # 1 is safe for concurrency
"""The maximum number of threads that PyTables should use internally in
Blosc.  If `None`, it is automatically set to the number of cores in
your machine.  For applications that use several PyTables instances
concurrently and so as to avoid locking problems, the recommended value
is 1.  In other cases a value of 2 or 4 could make sense.

"""

USER_BLOCK_SIZE = 0
"""Sets the user block size of a file.

The default user block size is 0; it may be set to any power of 2 equal
to 512 or greater (512, 1024, 2048, etc.).

.. versionadded:: 3.0

"""

ALLOW_PADDING = True
"""Allow padding in compound data types.

Starting on version 3.5 padding is honored during copies, or when tables
are created from NumPy structured arrays with padding (e.g. `align=True`).
If you actually want to get rid of any possible padding in new
datasets/attributes (i.e. the previous behaviour), set this to `False`.

.. versionadded:: 3.5

"""


# HDF5 driver management
# ----------------------
DRIVER = None
"""The HDF5 driver that should be used for reading/writing to the file.

Following drivers are supported:

    * H5FD_SEC2: this driver uses POSIX file-system functions like read
      and write to perform I/O to a single, permanent file on local
      disk with no system buffering.
      This driver is POSIX-compliant and is the default file driver for
      all systems.

    * H5FD_DIRECT: this is the H5FD_SEC2 driver except data is written
      to or read from the file synchronously without being cached by
      the system.

    * H5FD_WINDOWS: this driver was modified in HDF5-1.8.8 to be a
      wrapper of the POSIX driver, H5FD_SEC2. This change should not
      affect user applications.

    * H5FD_STDIO: this driver uses functions from the standard C
      stdio.h to perform I/O to a single, permanent file on local disk
      with additional system buffering.

    * H5FD_CORE: with this driver, an application can work with a file
      in memory for faster reads and writes. File contents are kept in
      memory until the file is closed. At closing, the memory version
      of the file can be written back to disk or abandoned.

    * H5FD_SPLIT: this file driver splits a file into two parts.
      One part stores metadata, and the other part stores raw data.
      This splitting a file into two parts is a limited case of the
      Multi driver.

The following drivers are not currently supported:

    * H5FD_LOG: this is the H5FD_SEC2 driver with logging capabilities.

    * H5FD_FAMILY: with this driver, the HDF5 file’s address space is
      partitioned into pieces and sent to separate storage files using
      an underlying driver of the user’s choice.
      This driver is for systems that do not support files larger than
      2 gigabytes.

    * H5FD_MULTI: with this driver, data can be stored in multiple
      files according to the type of the data. I/O might work better if
      data is stored in separate files based on the type of data.
      The Split driver is a special case of this driver.

    * H5FD_MPIO: this is the standard HDF5 file driver for parallel
      file systems. This driver uses the MPI standard for both
      communication and file I/O.

    * H5FD_MPIPOSIX: this parallel file system driver uses MPI for
      communication and POSIX file-system calls for file I/O.

    * H5FD_STREAM: this driver is no longer available.

.. seealso:: the `Drivers section`_ of the `HDF5 User's Guide`_ for
   more information.

.. note::

    not all supported drivers are always available. For example the
    H5FD_WINDOWS driver is not available on non Windows platforms.

    If the user try to use a driver that is not available on the target
    platform a :exc:`RuntimeError` is raised.

.. versionadded:: 3.0

.. _`Drivers section`:
    http://www.hdfgroup.org/HDF5/doc/UG/08_TheFile.html#Drivers
.. _`HDF5 User's Guide`: http://www.hdfgroup.org/HDF5/doc/UG/index.html

"""

DRIVER_DIRECT_ALIGNMENT = 0
"""Specifies the required alignment boundary in memory.

A value of 0 (zero) means to use HDF5 Library’s default value.

.. versionadded:: 3.0

"""

DRIVER_DIRECT_BLOCK_SIZE = 0
"""Specifies the file system block size.

A value of 0 (zero) means to use HDF5 Library’s default value of 4KB.

.. versionadded:: 3.0

"""

DRIVER_DIRECT_CBUF_SIZE = 0
"""Specifies the copy buffer size.

A value of 0 (zero) means to use HDF5 Library’s default value.

.. versionadded:: 3.0

"""

# DRIVER_LOG_FLAGS = 0x0001ffff
# """Flags specifying the types of logging activity.
#
# .. versionadded:: 3.0
#
# .. seeealso::
#     http://www.hdfgroup.org/HDF5/doc/RM/RM_H5P.html#Property-SetFaplLog
#
# """
#
#  DRIVER_LOG_BUF_SIZE = 4 * _KB
# """The size of the logging buffers, in bytes.
#
#  One buffer of size DRIVER_LOG_BUF_SIZE will be created for each of
#  H5FD_LOG_FILE_READ, H5FD_LOG_FILE_WRITE and H5FD_LOG_FLAVOR when those
#  flags are set; these buffers will not grow as the file increases in
#  size.
#
# .. versionadded:: 3.0
#
# """

DRIVER_CORE_INCREMENT = 64 * _KB
"""Core driver memory increment.

Specifies the increment by which allocated memory is to be increased
each time more memory is required.

.. versionadded:: 3.0

"""

DRIVER_CORE_BACKING_STORE = 1
"""Enables backing store for the core driver.

With the H5FD_CORE driver, if the DRIVER_CORE_BACKING_STORE is set
to 1 (True), the file contents are flushed to a file with the same name
as this core file when the file is closed or access to the file is
terminated in memory.

The application is allowed to open an existing file with H5FD_CORE
driver. In that case, if the DRIVER_CORE_BACKING_STORE is set to 1 and
the flags for :func:`tables.open_file` is set to H5F_ACC_RDWR, any change
to the file contents are saved to the file when the file is closed.
If backing_store is set to 0 and the flags for :func:`tables.open_file`
is set to H5F_ACC_RDWR, any change to the file contents will be lost
when the file is closed. If the flags for :func:`tables.open_file` is
set to H5F_ACC_RDONLY, no change to the file is allowed either in
memory or on file.

.. versionadded:: 3.0

"""

DRIVER_CORE_IMAGE = None
"""String containing an HDF5 file image.

If this option is passed to the :func:`tables.open_file` function then the
returned file object is set up using the specified image.

A file image can be retrieved from an existing (and opened) file object
using the :meth:`tables.File.get_file_image` method.

.. versionadded:: 3.0

"""

DRIVER_SPLIT_META_EXT = "-m.h5"
"""The extension for the metadata file used by the H5FD_SPLIT driver.

If this option is passed to the :func:`tables.openFile` function along
with driver='H5FD_SPLIT', the extension is appended to the name passed
as the first parameter to form the name of the metadata file. If the
string '%s' is used in the extension, the metadata file name is formed
by replacing '%s' with the name passed as the first parameter instead.

.. versionadded:: 3.1

"""

DRIVER_SPLIT_RAW_EXT = "-r.h5"
"""The extension for the raw data file used by the H5FD_SPLIT driver.

If this option is passed to the :func:`tables.openFile` function along
with driver='H5FD_SPLIT', the extension is appended to the name passed
as the first parameter to form the name of the raw data file. If the
string '%s' is used in the extension, the raw data file name is formed
by replacing '%s' with the name passed as the first parameter instead.

.. versionadded:: 3.1

"""
