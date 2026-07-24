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

"""
FileSystem abstraction to interact with various local and remote filesystems.
"""

from pyarrow.util import _is_path_like, _stringify_path

from pyarrow._fs import (  # noqa
    FileSelector,
    FileType,
    FileInfo,
    FileSystem,
    LocalFileSystem,
    SubTreeFileSystem,
    _MockFileSystem,
    FileSystemHandler,
    PyFileSystem,
    _copy_files,
    _copy_files_selector,
)

# For backward compatibility.
FileStats = FileInfo

_not_imported = []
try:
    from pyarrow._azurefs import AzureFileSystem  # noqa
except ImportError:
    _not_imported.append("AzureFileSystem")

try:
    from pyarrow._hdfs import HadoopFileSystem  # noqa
except ImportError:
    _not_imported.append("HadoopFileSystem")

try:
    from pyarrow._gcsfs import GcsFileSystem  # noqa
except ImportError:
    _not_imported.append("GcsFileSystem")

try:
    from pyarrow._s3fs import (  # noqa
        AwsDefaultS3RetryStrategy, AwsStandardS3RetryStrategy,
        S3FileSystem, S3LogLevel, S3RetryStrategy, ensure_s3_initialized,
        finalize_s3, ensure_s3_finalized, initialize_s3, resolve_s3_region)
except ImportError:
    _not_imported.append("S3FileSystem")
else:
    # GH-38364: we don't initialize S3 eagerly as that could lead
    # to crashes at shutdown even when S3 isn't used.
    # Instead, S3 is initialized lazily using `ensure_s3_initialized`
    # in assorted places.
    import atexit
    atexit.register(ensure_s3_finalized)


def __getattr__(name):
    if name in _not_imported:
        raise ImportError(
            "The pyarrow installation is not built with support for "
            f"'{name}'"
        )

    raise AttributeError(
        f"module 'pyarrow.fs' has no attribute '{name}'"
    )


def _ensure_filesystem(filesystem, *, use_mmap=False):
    if isinstance(filesystem, FileSystem):
        return filesystem
    elif isinstance(filesystem, str):
        # create a filesystem from a URI string, note that the `path` part of the URI
        # is treated as a prefix if specified, so the filesystem is wrapped in a
        # SubTreeFileSystem
        if use_mmap:
            raise ValueError(
                "Specifying to use memory mapping not supported for "
                "filesystem specified as an URI string"
            )
        fs, path = FileSystem.from_uri(filesystem)
        prefix = fs.normalize_path(path)
        if prefix:
            # validate that the prefix is pointing to a directory
            prefix_info = fs.get_file_info([prefix])[0]
            if prefix_info.type != FileType.Directory:
                raise ValueError(
                    "The path component of the filesystem URI must point to a "
                    f"directory but it has a type: `{prefix_info.type.name}`. The path "
                    f"component is `{prefix_info.path}` and the given filesystem URI "
                    f"is `{filesystem}`"
                )
            fs = SubTreeFileSystem(prefix, fs)
        return fs
    else:
        # handle fsspec-compatible filesystems
        try:
            import fsspec
        except ImportError:
            pass
        else:
            if isinstance(filesystem, fsspec.AbstractFileSystem):
                if type(filesystem).__name__ == 'LocalFileSystem':
                    # In case its a simple LocalFileSystem, use native arrow one
                    return LocalFileSystem(use_mmap=use_mmap)
                return PyFileSystem(FSSpecHandler(filesystem))

        raise TypeError(
            f"Unrecognized filesystem: {type(filesystem)}. `filesystem` argument must "
            "be a FileSystem instance or a valid file system URI"
        )


def _resolve_filesystem_and_path(path, filesystem=None, *, memory_map=False):
    """
    Return filesystem/path from path which could be an URI or a plain
    filesystem path or a combination of fsspec protocol and URI.
    """
    if not _is_path_like(path):
        if filesystem is not None:
            raise ValueError(
                "'filesystem' passed but the specified path is file-like, so"
                " there is nothing to open with 'filesystem'."
            )
        return filesystem, path

    if filesystem is not None:
        filesystem = _ensure_filesystem(filesystem, use_mmap=memory_map)
        if isinstance(filesystem, LocalFileSystem):
            path = _stringify_path(path)
        elif not isinstance(path, str):
            raise TypeError(
                "Expected string path; path-like objects are only allowed "
                "with a local filesystem"
            )
        path = filesystem.normalize_path(path)
        return filesystem, path

    path = _stringify_path(path)

    # if filesystem is not given, try to automatically determine one
    # first check if the file exists as a local (relative) file path
    # if not then try to parse the path as an URI
    filesystem = LocalFileSystem(use_mmap=memory_map)

    try:
        file_info = filesystem.get_file_info(path)
    except ValueError:  # ValueError means path is likely an URI
        file_info = None
        exists_locally = False
    else:
        exists_locally = (file_info.type != FileType.NotFound)

    # if the file or directory doesn't exists locally, then assume that
    # the path is an URI describing the file system as well
    if not exists_locally:
        try:
            filesystem, path = FileSystem.from_uri(path)
        except ValueError as e:
            msg = str(e)
            if "empty scheme" in msg or "Cannot parse URI" in msg:
                # neither an URI nor a locally existing path, so assume that
                # local path was given and propagate a nicer file not found
                # error instead of a more confusing scheme parsing error
                pass
            else:
                raise e
    else:
        path = filesystem.normalize_path(path)

    return filesystem, path


def copy_files(source, destination,
               source_filesystem=None, destination_filesystem=None,
               *, chunk_size=1024*1024, use_threads=True):
    """
    Copy files between FileSystems.

    This functions allows you to recursively copy directories of files from
    one file system to another, such as from S3 to your local machine.

    Parameters
    ----------
    source : string
        Source file path or URI to a single file or directory.
        If a directory, files will be copied recursively from this path.
    destination : string
        Destination file path or URI. If `source` is a file, `destination`
        is also interpreted as the destination file (not directory).
        Directories will be created as necessary.
    source_filesystem : FileSystem, optional
        Source filesystem, needs to be specified if `source` is not a URI,
        otherwise inferred.
    destination_filesystem : FileSystem, optional
        Destination filesystem, needs to be specified if `destination` is not
        a URI, otherwise inferred.
    chunk_size : int, default 1MB
        The maximum size of block to read before flushing to the
        destination file. A larger chunk_size will use more memory while
        copying but may help accommodate high latency FileSystems.
    use_threads : bool, default True
        Whether to use multiple threads to accelerate copying.

    Examples
    --------
    Inspect an S3 bucket's files:

    >>> s3, path = fs.FileSystem.from_uri(
    ...            "s3://registry.opendata.aws/roda/ndjson/")
    >>> selector = fs.FileSelector(path)
    >>> s3.get_file_info(selector)
    [<FileInfo for 'registry.opendata.aws/roda/ndjson/index.ndjson':...]

    Copy one file from S3 bucket to a local directory:

    >>> fs.copy_files("s3://registry.opendata.aws/roda/ndjson/index.ndjson",
    ...               f"file:///{local_path}/index_copy.ndjson")

    >>> fs.LocalFileSystem().get_file_info(str(local_path)+
    ...                                    '/index_copy.ndjson')
    <FileInfo for '.../index_copy.ndjson': type=FileType.File, size=...>

    Copy file using a FileSystem object:

    >>> fs.copy_files("registry.opendata.aws/roda/ndjson/index.ndjson",
    ...               f"file:///{local_path}/index_copy.ndjson",
    ...               source_filesystem=fs.S3FileSystem())
    """
    source_fs, source_path = _resolve_filesystem_and_path(
        source, source_filesystem
    )
    destination_fs, destination_path = _resolve_filesystem_and_path(
        destination, destination_filesystem
    )

    file_info = source_fs.get_file_info(source_path)
    if file_info.type == FileType.Directory:
        source_sel = FileSelector(source_path, recursive=True)
        _copy_files_selector(source_fs, source_sel,
                             destination_fs, destination_path,
                             chunk_size, use_threads)
    else:
        _copy_files(source_fs, source_path,
                    destination_fs, destination_path,
                    chunk_size, use_threads)


class FSSpecHandler(FileSystemHandler):
    """
    Handler for fsspec-based Python filesystems.

    https://filesystem-spec.readthedocs.io/en/latest/index.html

    Parameters
    ----------
    fs : FSSpec-compliant filesystem instance

    Examples
    --------
    >>> PyFileSystem(FSSpecHandler(fsspec_fs)) # doctest: +SKIP
    """

    def __init__(self, fs):
        self.fs = fs

    def __eq__(self, other):
        if isinstance(other, FSSpecHandler):
            return self.fs == other.fs
        return NotImplemented

    def __ne__(self, other):
        if isinstance(other, FSSpecHandler):
            return self.fs != other.fs
        return NotImplemented

    def get_type_name(self):
        protocol = self.fs.protocol
        if isinstance(protocol, list):
            protocol = protocol[0]
        return f"fsspec+{protocol}"

    def normalize_path(self, path):
        return path

    @staticmethod
    def _create_file_info(path, info):
        size = info["size"]
        if info["type"] == "file":
            ftype = FileType.File
        elif info["type"] == "directory":
            ftype = FileType.Directory
            # some fsspec filesystems include a file size for directories
            size = None
        else:
            ftype = FileType.Unknown
        return FileInfo(path, ftype, size=size, mtime=info.get("mtime", None))

    def get_file_info(self, paths):
        infos = []
        for path in paths:
            try:
                info = self.fs.info(path)
            except FileNotFoundError:
                infos.append(FileInfo(path, FileType.NotFound))
            else:
                infos.append(self._create_file_info(path, info))
        return infos

    def get_file_info_selector(self, selector):
        if not self.fs.isdir(selector.base_dir):
            if self.fs.exists(selector.base_dir):
                raise NotADirectoryError(selector.base_dir)
            else:
                if selector.allow_not_found:
                    return []
                else:
                    raise FileNotFoundError(selector.base_dir)

        if selector.recursive:
            maxdepth = None
        else:
            maxdepth = 1

        infos = []
        selected_files = self.fs.find(
            selector.base_dir, maxdepth=maxdepth, withdirs=True, detail=True
        )
        for path, info in selected_files.items():
            _path = path.strip("/")
            base_dir = selector.base_dir.strip("/")
            # Need to exclude base directory from selected files if present
            # (fsspec filesystems, see GH-37555)
            if _path != base_dir:
                infos.append(self._create_file_info(path, info))

        return infos

    def create_dir(self, path, recursive):
        # mkdir also raises FileNotFoundError when base directory is not found
        try:
            self.fs.mkdir(path, create_parents=recursive)
        except FileExistsError:
            pass

    def delete_dir(self, path):
        self.fs.rm(path, recursive=True)

    def _delete_dir_contents(self, path, missing_dir_ok):
        try:
            subpaths = self.fs.listdir(path, detail=False)
        except FileNotFoundError:
            if missing_dir_ok:
                return
            raise
        for subpath in subpaths:
            if self.fs.isdir(subpath):
                self.fs.rm(subpath, recursive=True)
            elif self.fs.isfile(subpath):
                self.fs.rm(subpath)

    def delete_dir_contents(self, path, missing_dir_ok):
        if path.strip("/") == "":
            raise ValueError(
                "delete_dir_contents called on path '", path, "'")
        self._delete_dir_contents(path, missing_dir_ok)

    def delete_root_dir_contents(self):
        self._delete_dir_contents("/", False)

    def delete_file(self, path):
        # fs.rm correctly raises IsADirectoryError when `path` is a directory
        # instead of a file and `recursive` is not set to True
        if not self.fs.exists(path):
            raise FileNotFoundError(path)
        self.fs.rm(path)

    def move(self, src, dest):
        self.fs.mv(src, dest, recursive=True)

    def copy_file(self, src, dest):
        # fs.copy correctly raises IsADirectoryError when `src` is a directory
        # instead of a file
        self.fs.copy(src, dest)

    # TODO can we read/pass metadata (e.g. Content-Type) in the methods below?

    def open_input_stream(self, path):
        from pyarrow import PythonFile

        if not self.fs.isfile(path):
            raise FileNotFoundError(path)

        return PythonFile(self.fs.open(path, mode="rb"), mode="r")

    def open_input_file(self, path):
        from pyarrow import PythonFile

        if not self.fs.isfile(path):
            raise FileNotFoundError(path)

        return PythonFile(self.fs.open(path, mode="rb"), mode="r")

    def open_output_stream(self, path, metadata):
        from pyarrow import PythonFile

        return PythonFile(self.fs.open(path, mode="wb"), mode="w")

    def open_append_stream(self, path, metadata):
        from pyarrow import PythonFile

        return PythonFile(self.fs.open(path, mode="ab"), mode="w")
