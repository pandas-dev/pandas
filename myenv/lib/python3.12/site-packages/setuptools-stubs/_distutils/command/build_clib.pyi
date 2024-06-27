from _typeshed import Incomplete

from ..cmd import Command

class build_clib(Command):
    description: str
    user_options: Incomplete
    boolean_options: Incomplete
    help_options: Incomplete
    build_clib: Incomplete
    build_temp: Incomplete
    libraries: Incomplete
    include_dirs: Incomplete
    define: Incomplete
    undef: Incomplete
    debug: Incomplete
    force: int
    compiler: Incomplete
    def initialize_options(self) -> None: ...
    def finalize_options(self) -> None: ...
    def run(self) -> None: ...
    def check_library_list(self, libraries) -> None: ...
    def get_library_names(self): ...
    def get_source_files(self): ...
    def build_libraries(self, libraries) -> None: ...
