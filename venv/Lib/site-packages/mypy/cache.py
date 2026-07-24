"""
This module contains high-level logic for fixed format serialization.

Lower-level parts are implemented in C in mypyc/lib-rt/internal/librt_internal.c
Short summary of low-level functionality:
* integers are automatically serialized as 1, 2, or 4 bytes, or arbitrary length.
* str/bytes are serialized as size (1, 2, or 4 bytes) followed by bytes buffer.
* floats are serialized as C doubles.

At high-level we add type tags as needed so that our format is self-descriptive.
More precisely:
* False, True, and None are stored as just a tag: 0, 1, 2 correspondingly.
* builtin primitives like int/str/bytes/float are stored as their type tag followed
  by bare (low-level) representation of the value. Reserved tag range for primitives is
  3 ... 19.
* generic (heterogeneous) list are stored as tag, followed by bare size, followed by
  sequence of tagged values.
* homogeneous lists of primitives are stored as tag, followed by bare size, followed
  by sequence of bare values.
* reserved tag range for sequence-like builtins is 20 ... 29
* currently we have only one mapping-like format: string-keyed dictionary with heterogeneous
  values. It is stored as tag, followed by bare size, followed by sequence of pairs: bare
  string key followed by tagged value.
* reserved tag range for mapping-like builtins is 30 ... 39
* there is an additional reserved tag range 40 ... 49 for any other builtin collections.
* custom classes (like types, symbols etc.) are stored as tag, followed by a sequence of
  tagged field values, followed by a special end tag 255. Names of class fields are
  *not* stored, the caller should know the field names and order for the given class tag.
* reserved tag range for symbols (TypeInfo, Var, etc) is 50 ... 79.
* class Instance is the only exception from the above format (since it is the most common one).
  It has two extra formats: few most common instances like "builtins.object" are stored as
  instance tag followed by a secondary tag, other plain non-generic instances are stored as
  instance tag followed by secondary tag followed by fullname as bare string. All generic
  readers must handle these.
* reserved tag range for Instance type formats is 80 ... 99, for other types it is 100 ... 149.
* tag 254 is reserved for if we would ever need to extend the tag range to indicated second tag
  page. Tags 150 ... 253 are free for everything else (e.g. AST nodes etc).

General convention is that custom classes implement write() and read() methods for FF
serialization. The write method should write both class tag and end tag. The read method
conventionally *does not* read the start tag (to simplify logic for unions). Known exceptions
are MypyFile.read() and SymbolTableNode.read(), since those two never appear in a union.

If any of these details change, or if the structure of CacheMeta changes please
bump CACHE_VERSION below.
"""

from __future__ import annotations

from collections.abc import Sequence
from typing import Any, Final, TypeAlias as _TypeAlias

from librt.internal import (
    ReadBuffer as ReadBuffer,
    WriteBuffer as WriteBuffer,
    read_bool as read_bool,
    read_bytes as read_bytes_bare,
    read_float as read_float_bare,
    read_int as read_int_bare,
    read_str as read_str_bare,
    read_tag as read_tag,
    write_bool as write_bool,
    write_bytes as write_bytes_bare,
    write_float as write_float_bare,
    write_int as write_int_bare,
    write_str as write_str_bare,
    write_tag as write_tag,
)
from mypy_extensions import u8

# High-level cache layout format
CACHE_VERSION: Final = 8

# Type used internally to represent errors:
#   (path, line, column, end_line, end_column, severity, message, code)
ErrorTuple: _TypeAlias = tuple[str | None, int, int, int, int, str, str, str | None]


class CacheMeta:
    """Class representing cache metadata for a module.

    This class represents the data known after checking module interface only, i.e.
    this doesn't have: error messages and indirect dependencies, these are stored
    in CacheMetaEx.
    """

    def __init__(
        self,
        *,
        id: str,
        path: str,
        mtime: int,
        size: int,
        hash: str,
        dependencies: list[str],
        data_mtime: int,
        data_file: str,
        suppressed: list[str],
        imports_ignored: dict[int, list[str]],
        options: dict[str, object],
        suppressed_deps_opts: bytes,
        dep_prios: list[int],
        dep_lines: list[int],
        dep_hashes: list[bytes],
        interface_hash: bytes,
        trans_dep_hash: bytes,
        version_id: str,
        ignore_all: bool,
        plugin_data: Any,
    ) -> None:
        self.id = id
        self.path = path
        self.mtime = mtime  # source file mtime
        self.size = size  # source file size
        self.hash = hash  # source file hash (as a hex string for historical reasons)
        self.dependencies = dependencies  # names of imported modules
        self.data_mtime = data_mtime  # mtime of data_file
        self.data_file = data_file  # path of <id>.data.json or <id>.data.ff
        self.suppressed = suppressed  # dependencies that weren't imported
        self.imports_ignored = imports_ignored  # type ignore codes by line
        self.options = options  # build options snapshot
        self.suppressed_deps_opts = suppressed_deps_opts  # hash of import-related options
        # dep_prios and dep_lines are both aligned with dependencies + suppressed
        self.dep_prios = dep_prios
        self.dep_lines = dep_lines
        # dep_hashes list is aligned with dependencies only
        self.dep_hashes = dep_hashes  # list of interface_hash for dependencies
        self.interface_hash = interface_hash  # hash representing the public interface
        self.trans_dep_hash = trans_dep_hash  # hash of import structure (transitive)
        self.version_id = version_id  # mypy version for cache invalidation
        self.ignore_all = ignore_all  # if errors were ignored
        self.plugin_data = plugin_data  # config data from plugins

    def serialize(self) -> dict[str, Any]:
        return {
            "id": self.id,
            "path": self.path,
            "mtime": self.mtime,
            "size": self.size,
            "hash": self.hash,
            "data_mtime": self.data_mtime,
            "dependencies": self.dependencies,
            "suppressed": self.suppressed,
            "imports_ignored": {str(line): codes for line, codes in self.imports_ignored.items()},
            "options": self.options,
            "suppressed_deps_opts": self.suppressed_deps_opts.hex(),
            "dep_prios": self.dep_prios,
            "dep_lines": self.dep_lines,
            "dep_hashes": [dep.hex() for dep in self.dep_hashes],
            "interface_hash": self.interface_hash.hex(),
            "trans_dep_hash": self.trans_dep_hash.hex(),
            "version_id": self.version_id,
            "ignore_all": self.ignore_all,
            "plugin_data": self.plugin_data,
        }

    @classmethod
    def deserialize(cls, meta: dict[str, Any], data_file: str) -> CacheMeta | None:
        try:
            return CacheMeta(
                id=meta["id"],
                path=meta["path"],
                mtime=meta["mtime"],
                size=meta["size"],
                hash=meta["hash"],
                dependencies=meta["dependencies"],
                data_mtime=meta["data_mtime"],
                data_file=data_file,
                suppressed=meta["suppressed"],
                imports_ignored={
                    int(line): codes for line, codes in meta["imports_ignored"].items()
                },
                options=meta["options"],
                suppressed_deps_opts=bytes.fromhex(meta["suppressed_deps_opts"]),
                dep_prios=meta["dep_prios"],
                dep_lines=meta["dep_lines"],
                dep_hashes=[bytes.fromhex(dep) for dep in meta["dep_hashes"]],
                interface_hash=bytes.fromhex(meta["interface_hash"]),
                trans_dep_hash=bytes.fromhex(meta["trans_dep_hash"]),
                version_id=meta["version_id"],
                ignore_all=meta["ignore_all"],
                plugin_data=meta["plugin_data"],
            )
        except (KeyError, ValueError):
            return None

    def write(self, data: WriteBuffer) -> None:
        write_str(data, self.id)
        write_str(data, self.path)
        write_int(data, self.mtime)
        write_int(data, self.size)
        write_str(data, self.hash)
        write_str_list(data, self.dependencies)
        write_int(data, self.data_mtime)
        write_str_list(data, self.suppressed)
        write_int_bare(data, len(self.imports_ignored))
        for line, codes in self.imports_ignored.items():
            write_int(data, line)
            write_str_list(data, codes)
        write_json(data, self.options)
        write_bytes(data, self.suppressed_deps_opts)
        write_int_list(data, self.dep_prios)
        write_int_list(data, self.dep_lines)
        write_bytes_list(data, self.dep_hashes)
        write_bytes(data, self.interface_hash)
        write_bytes(data, self.trans_dep_hash)
        write_str(data, self.version_id)
        write_bool(data, self.ignore_all)
        # Plugin data may be not a dictionary, so we use
        # a more generic write_json_value() here.
        write_json_value(data, self.plugin_data)

    @classmethod
    def read(cls, data: ReadBuffer, data_file: str) -> CacheMeta | None:
        try:
            return CacheMeta(
                id=read_str(data),
                path=read_str(data),
                mtime=read_int(data),
                size=read_int(data),
                hash=read_str(data),
                dependencies=read_str_list(data),
                data_mtime=read_int(data),
                data_file=data_file,
                suppressed=read_str_list(data),
                imports_ignored={
                    read_int(data): read_str_list(data) for _ in range(read_int_bare(data))
                },
                options=read_json(data),
                suppressed_deps_opts=read_bytes(data),
                dep_prios=read_int_list(data),
                dep_lines=read_int_list(data),
                dep_hashes=read_bytes_list(data),
                interface_hash=read_bytes(data),
                trans_dep_hash=read_bytes(data),
                version_id=read_str(data),
                ignore_all=read_bool(data),
                plugin_data=read_json_value(data),
            )
        except (ValueError, AssertionError):
            return None


class CacheMetaEx:
    """Class representing "implementation-specific" part of cache metadata for a module."""

    def __init__(
        self,
        dependencies: list[str],
        suppressed: list[str],
        dep_hashes: list[bytes],
        error_lines: list[ErrorTuple],
    ) -> None:
        self.dependencies = dependencies
        self.suppressed = suppressed
        self.dep_hashes = dep_hashes
        self.error_lines = error_lines

    def serialize(self) -> dict[str, Any]:
        return {
            "dependencies": self.dependencies,
            "suppressed": self.suppressed,
            "dep_hashes": [dep.hex() for dep in self.dep_hashes],
            "error_lines": self.error_lines,
        }

    @classmethod
    def deserialize(cls, meta: dict[str, Any]) -> CacheMetaEx | None:
        try:
            return CacheMetaEx(
                dependencies=meta["dependencies"],
                suppressed=meta["suppressed"],
                dep_hashes=[bytes.fromhex(dep) for dep in meta["dep_hashes"]],
                error_lines=[tuple(err) for err in meta["error_lines"]],
            )
        except (KeyError, ValueError):
            return None

    def write(self, data: WriteBuffer) -> None:
        write_str_list(data, self.dependencies)
        write_str_list(data, self.suppressed)
        write_bytes_list(data, self.dep_hashes)
        write_errors(data, self.error_lines)

    @classmethod
    def read(cls, data: ReadBuffer) -> CacheMetaEx | None:
        try:
            return CacheMetaEx(
                dependencies=read_str_list(data),
                suppressed=read_str_list(data),
                dep_hashes=read_bytes_list(data),
                error_lines=read_errors(data),
            )
        except (ValueError, AssertionError):
            return None


# Always use this type alias to refer to type tags.
Tag = u8

# Note: all tags should be kept in sync with lib-rt/internal/librt_internal.c.
# Primitives.
LITERAL_FALSE: Final[Tag] = 0
LITERAL_TRUE: Final[Tag] = 1
LITERAL_NONE: Final[Tag] = 2
LITERAL_INT: Final[Tag] = 3
LITERAL_STR: Final[Tag] = 4
LITERAL_BYTES: Final[Tag] = 5
LITERAL_FLOAT: Final[Tag] = 6
LITERAL_COMPLEX: Final[Tag] = 7

# Collections.
LIST_GEN: Final[Tag] = 20
LIST_INT: Final[Tag] = 21
LIST_STR: Final[Tag] = 22
LIST_BYTES: Final[Tag] = 23
TUPLE_GEN: Final[Tag] = 24
DICT_STR_GEN: Final[Tag] = 30
DICT_INT_GEN: Final[Tag] = 31

# Misc classes.
EXTRA_ATTRS: Final[Tag] = 150
DT_SPEC: Final[Tag] = 151
# Four integers representing source file (line, column) range.
LOCATION: Final[Tag] = 152

RESERVED: Final[Tag] = 254
END_TAG: Final[Tag] = 255


def read_literal(data: ReadBuffer, tag: Tag) -> int | str | bool | float:
    if tag == LITERAL_INT:
        return read_int_bare(data)
    elif tag == LITERAL_STR:
        return read_str_bare(data)
    elif tag == LITERAL_FALSE:
        return False
    elif tag == LITERAL_TRUE:
        return True
    elif tag == LITERAL_FLOAT:
        return read_float_bare(data)
    assert False, f"Unknown literal tag {tag}"


# There is an intentional asymmetry between read and write for literals because
# None and/or complex values are only allowed in some contexts but not in others.
def write_literal(data: WriteBuffer, value: int | str | bool | float | complex | None) -> None:
    if isinstance(value, bool):
        write_bool(data, value)
    elif isinstance(value, int):
        write_tag(data, LITERAL_INT)
        write_int_bare(data, value)
    elif isinstance(value, str):
        write_tag(data, LITERAL_STR)
        write_str_bare(data, value)
    elif isinstance(value, float):
        write_tag(data, LITERAL_FLOAT)
        write_float_bare(data, value)
    elif isinstance(value, complex):
        write_tag(data, LITERAL_COMPLEX)
        write_float_bare(data, value.real)
        write_float_bare(data, value.imag)
    else:
        write_tag(data, LITERAL_NONE)


def read_int(data: ReadBuffer) -> int:
    assert read_tag(data) == LITERAL_INT
    return read_int_bare(data)


def write_int(data: WriteBuffer, value: int) -> None:
    write_tag(data, LITERAL_INT)
    write_int_bare(data, value)


def read_str(data: ReadBuffer) -> str:
    assert read_tag(data) == LITERAL_STR
    return read_str_bare(data)


def write_str(data: WriteBuffer, value: str) -> None:
    write_tag(data, LITERAL_STR)
    write_str_bare(data, value)


def read_bytes(data: ReadBuffer) -> bytes:
    assert read_tag(data) == LITERAL_BYTES
    return read_bytes_bare(data)


def write_bytes(data: WriteBuffer, value: bytes) -> None:
    write_tag(data, LITERAL_BYTES)
    write_bytes_bare(data, value)


def read_int_opt(data: ReadBuffer) -> int | None:
    tag = read_tag(data)
    if tag == LITERAL_NONE:
        return None
    assert tag == LITERAL_INT
    return read_int_bare(data)


def write_int_opt(data: WriteBuffer, value: int | None) -> None:
    if value is not None:
        write_tag(data, LITERAL_INT)
        write_int_bare(data, value)
    else:
        write_tag(data, LITERAL_NONE)


def read_str_opt(data: ReadBuffer) -> str | None:
    tag = read_tag(data)
    if tag == LITERAL_NONE:
        return None
    assert tag == LITERAL_STR
    return read_str_bare(data)


def write_str_opt(data: WriteBuffer, value: str | None) -> None:
    if value is not None:
        write_tag(data, LITERAL_STR)
        write_str_bare(data, value)
    else:
        write_tag(data, LITERAL_NONE)


def read_int_list(data: ReadBuffer) -> list[int]:
    assert read_tag(data) == LIST_INT
    size = read_int_bare(data)
    return [read_int_bare(data) for _ in range(size)]


def write_int_list(data: WriteBuffer, value: list[int]) -> None:
    write_tag(data, LIST_INT)
    write_int_bare(data, len(value))
    for item in value:
        write_int_bare(data, item)


def read_str_list(data: ReadBuffer) -> list[str]:
    assert read_tag(data) == LIST_STR
    size = read_int_bare(data)
    return [read_str_bare(data) for _ in range(size)]


def write_str_list(data: WriteBuffer, value: Sequence[str]) -> None:
    write_tag(data, LIST_STR)
    write_int_bare(data, len(value))
    for item in value:
        write_str_bare(data, item)


def read_bytes_list(data: ReadBuffer) -> list[bytes]:
    assert read_tag(data) == LIST_BYTES
    size = read_int_bare(data)
    return [read_bytes_bare(data) for _ in range(size)]


def write_bytes_list(data: WriteBuffer, value: Sequence[bytes]) -> None:
    write_tag(data, LIST_BYTES)
    write_int_bare(data, len(value))
    for item in value:
        write_bytes_bare(data, item)


def read_str_opt_list(data: ReadBuffer) -> list[str | None]:
    assert read_tag(data) == LIST_GEN
    size = read_int_bare(data)
    return [read_str_opt(data) for _ in range(size)]


def write_str_opt_list(data: WriteBuffer, value: list[str | None]) -> None:
    write_tag(data, LIST_GEN)
    write_int_bare(data, len(value))
    for item in value:
        write_str_opt(data, item)


Value: _TypeAlias = None | int | str | bool

# Our JSON format is somewhat non-standard as we distinguish lists and tuples.
# This is convenient for some internal things, like mypyc plugin and error serialization.
JsonValue: _TypeAlias = (
    Value | list["JsonValue"] | dict[str, "JsonValue"] | tuple["JsonValue", ...]
)


def read_json_value(data: ReadBuffer) -> JsonValue:
    tag = read_tag(data)
    if tag == LITERAL_NONE:
        return None
    if tag == LITERAL_FALSE:
        return False
    if tag == LITERAL_TRUE:
        return True
    if tag == LITERAL_INT:
        return read_int_bare(data)
    if tag == LITERAL_STR:
        return read_str_bare(data)
    if tag == LIST_GEN:
        size = read_int_bare(data)
        return [read_json_value(data) for _ in range(size)]
    if tag == TUPLE_GEN:
        size = read_int_bare(data)
        return tuple(read_json_value(data) for _ in range(size))
    if tag == DICT_STR_GEN:
        size = read_int_bare(data)
        return {read_str_bare(data): read_json_value(data) for _ in range(size)}
    assert False, f"Invalid JSON tag: {tag}"


def write_json_value(data: WriteBuffer, value: JsonValue) -> None:
    if value is None:
        write_tag(data, LITERAL_NONE)
    elif isinstance(value, bool):
        write_bool(data, value)
    elif isinstance(value, int):
        write_tag(data, LITERAL_INT)
        write_int_bare(data, value)
    elif isinstance(value, str):
        write_tag(data, LITERAL_STR)
        write_str_bare(data, value)
    elif isinstance(value, list):
        write_tag(data, LIST_GEN)
        write_int_bare(data, len(value))
        for val in value:
            write_json_value(data, val)
    elif isinstance(value, tuple):
        write_tag(data, TUPLE_GEN)
        write_int_bare(data, len(value))
        for val in value:
            write_json_value(data, val)
    elif isinstance(value, dict):
        write_tag(data, DICT_STR_GEN)
        write_int_bare(data, len(value))
        for key in sorted(value):
            write_str_bare(data, key)
            write_json_value(data, value[key])
    else:
        assert False, f"Invalid JSON value: {value}"


# These are functions for JSON *dictionaries* specifically. Unfortunately, we
# must use imprecise types here, because the callers use imprecise types.
def read_json(data: ReadBuffer) -> dict[str, Any]:
    assert read_tag(data) == DICT_STR_GEN
    size = read_int_bare(data)
    return {read_str_bare(data): read_json_value(data) for _ in range(size)}


def write_json(data: WriteBuffer, value: dict[str, Any]) -> None:
    write_tag(data, DICT_STR_GEN)
    write_int_bare(data, len(value))
    for key in sorted(value):
        write_str_bare(data, key)
        write_json_value(data, value[key])


def write_errors(data: WriteBuffer, errs: list[ErrorTuple]) -> None:
    write_tag(data, LIST_GEN)
    write_int_bare(data, len(errs))
    for path, line, column, end_line, end_column, severity, message, code in errs:
        write_tag(data, TUPLE_GEN)
        write_str_opt(data, path)
        write_int(data, line)
        write_int(data, column)
        write_int(data, end_line)
        write_int(data, end_column)
        write_str(data, severity)
        write_str(data, message)
        write_str_opt(data, code)


def read_errors(data: ReadBuffer) -> list[ErrorTuple]:
    assert read_tag(data) == LIST_GEN
    result = []
    for _ in range(read_int_bare(data)):
        assert read_tag(data) == TUPLE_GEN
        result.append(
            (
                read_str_opt(data),
                read_int(data),
                read_int(data),
                read_int(data),
                read_int(data),
                read_str(data),
                read_str(data),
                read_str_opt(data),
            )
        )
    return result
