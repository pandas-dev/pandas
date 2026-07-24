from __future__ import annotations

from typing import Final


class Capsule:
    """Defines a C extension capsule that a primitive may require."""

    def __init__(self, name: str) -> None:
        # Module fullname, e.g. 'librt.base64'
        self.name: Final = name

    def __repr__(self) -> str:
        return f"Capsule(name={self.name!r})"

    def __eq__(self, other: object) -> bool:
        return isinstance(other, Capsule) and self.name == other.name

    def __hash__(self) -> int:
        return hash(("Capsule", self.name))

    def internal_dep(self) -> SourceDep:
        """Internal source dependency of the capsule that should only be included in the C extensions
        that depend on the capsule, eg. by importing a type or function from the capsule.
        """
        module = self.name.split(".")[-1]
        return SourceDep(f"{module}/librt_{module}_api.c", include_dirs=[module])

    def external_dep(self) -> HeaderDep:
        """External header dependency of the capsule that may be included in external headers of C
        extensions that depend on the capsule.

        The external headers of the C extensions are included by other C extensions that don't
        necessarily import the capsule. However, they may need type definitions from the capsule
        for types that are used in the exports table of the included C extensions.

        Only the external header should be included in this case because if the other C extension
        doesn't import the capsule, it also doesn't include the definition for its API table and
        including the internal header would result in undefined symbols.
        """
        module = self.name.split(".")[-1]
        return HeaderDep(f"{module}/librt_{module}.h", include_dirs=[module], internal=False)


class SourceDep:
    """Defines a C source file that a primitive may require.

    Each source file must also have a corresponding .h file (replace .c with .h)
    that gets implicitly #included if the source is used.
    include_dirs are passed to the C compiler when the file is compiled as a
    shared library separate from the C extension.
    """

    def __init__(
        self, path: str, *, include_dirs: list[str] | None = None, internal: bool = True
    ) -> None:
        # Relative path from mypyc/lib-rt, e.g. 'bytes_extra_ops.c'
        self.path: Final = path
        self.include_dirs: Final = include_dirs or []
        self.internal: Final = internal

    def __repr__(self) -> str:
        return f"SourceDep(path={self.path!r})"

    def __eq__(self, other: object) -> bool:
        return isinstance(other, SourceDep) and self.path == other.path

    def __hash__(self) -> int:
        return hash(("SourceDep", self.path))

    def get_header(self) -> str:
        """Get the header file path by replacing .c with .h"""
        return self.path.replace(".c", ".h")


class HeaderDep:
    """Defines a C header file that a primitive may require.

    The header gets explicitly #included if the dependency is used.
    include_dirs are passed to the C compiler when the generated extension
    is compiled separately and needs to include the header.
    """

    def __init__(
        self, path: str, *, include_dirs: list[str] | None = None, internal: bool = True
    ) -> None:
        # Relative path from mypyc/lib-rt, e.g. 'strings/librt_strings.h'
        self.path: Final = path
        self.include_dirs: Final = include_dirs or []
        self.internal: Final = internal

    def __repr__(self) -> str:
        return f"HeaderDep(path={self.path!r})"

    def __eq__(self, other: object) -> bool:
        return isinstance(other, HeaderDep) and self.path == other.path

    def __hash__(self) -> int:
        return hash(("HeaderDep", self.path))

    def get_header(self) -> str:
        return self.path


Dependency = Capsule | SourceDep | HeaderDep


LIBRT_STRINGS: Final = Capsule("librt.strings")
LIBRT_BASE64: Final = Capsule("librt.base64")
LIBRT_VECS: Final = Capsule("librt.vecs")
LIBRT_TIME: Final = Capsule("librt.time")
LIBRT_RANDOM: Final = Capsule("librt.random")

BYTES_EXTRA_OPS: Final = SourceDep("bytes_extra_ops.c")
BYTES_WRITER_EXTRA_OPS: Final = SourceDep("byteswriter_extra_ops.c")
STRING_WRITER_EXTRA_OPS: Final = SourceDep("stringwriter_extra_ops.c")
BYTEARRAY_EXTRA_OPS: Final = SourceDep("bytearray_extra_ops.c")
STR_EXTRA_OPS: Final = SourceDep("str_extra_ops.c")
VECS_EXTRA_OPS: Final = SourceDep("vecs_extra_ops.c")
