from __future__ import annotations

import sys


class CompilerOptions:
    def __init__(
        self,
        strip_asserts: bool = False,
        multi_file: bool = False,
        verbose: bool = False,
        separate: bool = False,
        target_dir: str | None = None,
        include_runtime_files: bool | None = None,
        capi_version: tuple[int, int] | None = None,
        python_version: tuple[int, int] | None = None,
        strict_dunder_typing: bool = False,
        group_name: str | None = None,
    ) -> None:
        self.strip_asserts = strip_asserts
        self.multi_file = multi_file
        self.verbose = verbose
        self.separate = separate
        self.global_opts = not separate
        self.target_dir = target_dir or "build"
        self.include_runtime_files = (
            include_runtime_files if include_runtime_files is not None else not multi_file
        )
        # The target Python C API version. Overriding this is mostly
        # useful in IR tests, since there's no guarantee that
        # binaries are backward compatible even if no recent API
        # features are used.
        self.capi_version = capi_version or sys.version_info[:2]
        self.python_version = python_version
        # Make possible to inline dunder methods in the generated code.
        # Typically, the convention is the dunder methods can return `NotImplemented`
        # even when its return type is just `bool`.
        # By enabling this option, this convention is no longer valid and the dunder
        # will assume the return type of the method strictly, which can lead to
        # more optimization opportunities.
        self.strict_dunders_typing = strict_dunder_typing
        # Override the automatic group name derived from the hash of module names.
        # This affects the names of generated .c, .h and shared library files.
        # This is only supported when compiling exactly one group, and a shared
        # library is generated (with shims). This can be used to make the output
        # file names more predictable.
        self.group_name = group_name
