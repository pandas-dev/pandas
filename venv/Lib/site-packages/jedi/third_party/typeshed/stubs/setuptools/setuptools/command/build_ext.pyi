from _typeshed import Incomplete
from typing import ClassVar

from setuptools.dist import Distribution

from .._distutils.command.build_ext import build_ext as _build_ext

have_rtld: bool
use_stubs: bool
libtype: str

def get_abi3_suffix(): ...

class build_ext(_build_ext):
    distribution: Distribution  # override distutils.dist.Distribution with setuptools.dist.Distribution
    editable_mode: ClassVar[bool]
    inplace: bool
    def run(self) -> None: ...
    def copy_extensions_to_source(self) -> None: ...
    def get_ext_filename(self, fullname): ...
    shlib_compiler: Incomplete
    shlibs: list[Incomplete]
    ext_map: dict[Incomplete, Incomplete]
    def initialize_options(self) -> None: ...
    extensions: list[Incomplete]
    def finalize_options(self) -> None: ...
    def setup_shlib_compiler(self) -> None: ...
    def get_export_symbols(self, ext): ...
    compiler: Incomplete
    def build_extension(self, ext) -> None: ...
    def links_to_dynamic(self, ext): ...
    def get_source_files(self) -> list[str]: ...
    def get_outputs(self) -> list[str]: ...
    def get_output_mapping(self) -> dict[str, str]: ...
    def write_stub(self, output_dir, ext, compile: bool = False) -> None: ...

def link_shared_object(
    self,
    objects,
    output_libname,
    output_dir=None,
    libraries=None,
    library_dirs=None,
    runtime_library_dirs=None,
    export_symbols=None,
    debug: bool = False,
    extra_preargs=None,
    extra_postargs=None,
    build_temp=None,
    target_lang=None,
) -> None: ...
