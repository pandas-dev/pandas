from typing import List, Optional, Tuple

class Extension:
    def __init__(
        self,
        name: str,
        sources: List[str],
        include_dirs: Optional[List[str]] = ...,
        define_macros: Optional[List[Tuple[str, Optional[str]]]] = ...,
        undef_macros: Optional[List[str]] = ...,
        library_dirs: Optional[List[str]] = ...,
        libraries: Optional[List[str]] = ...,
        runtime_library_dirs: Optional[List[str]] = ...,
        extra_objects: Optional[List[str]] = ...,
        extra_compile_args: Optional[List[str]] = ...,
        extra_link_args: Optional[List[str]] = ...,
        export_symbols: Optional[List[str]] = ...,
        swig_opts: Optional[str] = ...,  # undocumented
        depends: Optional[List[str]] = ...,
        language: Optional[str] = ...,
        optional: Optional[bool] = ...,
    ) -> None: ...
