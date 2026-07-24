from abc import abstractmethod
from collections.abc import Callable, Iterable
from threading import local
from types import ModuleType
from typing import Literal, NewType
from typing_extensions import TypeAlias

from fanstatic.compiler import Compiler, Minifier

_Renderer: TypeAlias = Callable[[str], str]

DEFAULT_SIGNATURE: str
VERSION_PREFIX: str
BUNDLE_PREFIX: str
NEEDED: str
DEBUG: str
MINIFIED: str

def set_resource_file_existence_checking(v: bool) -> None: ...
def set_auto_register_library(v: bool) -> None: ...

class UnknownResourceExtensionError(Exception): ...
class ModeResourceDependencyError(Exception): ...

UnknownResourceExtension = UnknownResourceExtensionError

class UnknownResourceError(Exception): ...
class ConfigurationError(Exception): ...
class LibraryDependencyCycleError(Exception): ...
class SlotError(Exception): ...

class Library:
    path: str
    name: str
    rootpath: str
    ignores: list[str]
    version: str | None
    known_resources: dict[str, Resource]
    known_assets: dict[str, Asset]
    module: ModuleType
    compilers: dict[str, Compiler]
    minifiers: dict[str, Minifier]
    def __init__(
        self,
        name: str,
        rootpath: str,
        ignores: list[str] | None = None,
        version: str | None = None,
        compilers: dict[str, Compiler] | None = None,
        minifiers: dict[str, Minifier] | None = None,
    ) -> None: ...
    def check_dependency_cycle(self, resource: Resource) -> None: ...
    def register(self, resource: Resource) -> None: ...
    def signature(self, recompute_hashes: bool = False, version_method: Callable[[str], str] | None = None) -> str: ...

def caller_dir() -> str: ...

class InclusionRenderers(dict[str, tuple[int, _Renderer]]):
    def register(self, extension: str, renderer: _Renderer, order: int | None = None) -> None: ...

inclusion_renderers: InclusionRenderers

def register_inclusion_renderer(extension: str, renderer: _Renderer, order: int | None = None) -> None: ...
def render_ico(url: str) -> str: ...
def render_css(url: str) -> str: ...
def render_js(url: str) -> str: ...
def render_print_css(url: str) -> str: ...
def render_screen_css(url: str) -> str: ...

class Renderable:
    @abstractmethod
    def render(self, library_url: str) -> str: ...

class Dependable:
    @property
    @abstractmethod
    def resources(self) -> set[Dependable]: ...
    @property
    @abstractmethod
    def depends(self) -> set[Dependable]: ...
    @property
    @abstractmethod
    def supports(self) -> set[Dependable]: ...
    def add_dependency(self, dependency: Dependable) -> None: ...
    @abstractmethod
    def set_dependencies(self, dependencies: Iterable[Dependable] | None) -> None: ...
    @abstractmethod
    def list_assets(self) -> set[Asset]: ...
    def list_supporting(self) -> set[Dependable]: ...

class Asset(Dependable):
    resources: set[Dependable]
    depends: set[Dependable]
    supports: set[Dependable]
    library: Library
    def __init__(self, library: Library, depends: Iterable[Dependable] | None = None) -> None: ...
    def set_dependencies(self, depends: Iterable[Dependable] | None) -> None: ...
    def list_assets(self) -> set[Asset]: ...

_NothingType = NewType("_NothingType", object)
NOTHING: _NothingType

class Resource(Renderable, Asset):
    relpath: str
    ext: str
    mode_parent: str | None
    compiler: Compiler
    source: str | None
    minifier: Minifier
    minified: Resource | None
    bottom: bool
    dont_bundle: bool
    renderer: _Renderer
    modes: dict[str, Resource]
    supersedes: list[Resource]
    rollups: list[Resource]
    def __init__(
        self,
        library: Library,
        relpath: str,
        depends: Iterable[Dependable] | None = None,
        supersedes: list[Resource] | None = None,
        bottom: bool = False,
        renderer: _Renderer | None = None,
        debug: str | Resource | None = None,
        dont_bundle: bool = False,
        minified: str | Resource | None = None,
        minifier: Minifier | _NothingType = ...,
        compiler: Compiler | _NothingType = ...,
        source: str | None = None,
        mode_parent: str | None = None,
    ) -> None: ...
    def fullpath(self, path: str | None = None) -> str: ...
    def compile(self, force: bool = False) -> None: ...
    def render(self, library_url: str) -> str: ...
    def mode(self, mode: str | None) -> Resource: ...
    def need(self, slots: dict[Slot, Resource] | None = None) -> None: ...

_RequiredDefaultMarkerType = NewType("_RequiredDefaultMarkerType", object)
REQUIRED_DEFAULT_MARKER: _RequiredDefaultMarkerType

class Slot(Asset):
    default: Resource | None
    ext: str
    required: bool
    def __init__(
        self,
        library: Library,
        extension: str,
        depends: Iterable[Dependable] | None = None,
        required: bool | _RequiredDefaultMarkerType = ...,
        default: Resource | None = None,
    ) -> None: ...

class FilledSlot(Renderable):
    filledby: Resource
    library: Library
    relpath: str
    bottom: bool
    rollups: list[Resource]
    dont_bundle: bool
    ext: str
    order: int
    renderer: _Renderer
    dependency_nr: int
    modes: dict[str, FilledSlot]
    def __init__(self, slot: Slot, resource: Resource) -> None: ...
    def render(self, library_url: str) -> str: ...
    def compile(self, force: bool = False) -> None: ...
    def mode(self, mode: str | None) -> FilledSlot: ...

class Group(Dependable):
    resources: set[Dependable]
    depends: set[Dependable]
    supports: set[Dependable]
    def __init__(self, depends: Iterable[Dependable]) -> None: ...
    def set_dependencies(self, depends: Iterable[Dependable]) -> None: ...  # type: ignore[override]
    def list_assets(self) -> set[Asset]: ...
    def need(self, slots: dict[Slot, Resource] | None = None) -> None: ...

GroupResource = Group

class NeededResources:
    def __init__(
        self,
        versioning: bool = False,
        versioning_use_md5: bool = False,
        recompute_hashes: bool = True,
        base_url: str | None = None,
        script_name: str | None = None,
        publisher_signature: str = "fanstatic",
        resources: Iterable[Dependable] | None = None,
    ) -> None: ...
    def has_resources(self) -> bool: ...
    def has_base_url(self) -> bool: ...
    def set_base_url(self, url: str) -> None: ...
    def need(self, resource: Resource | Group, slots: dict[Slot, Resource] | None = None) -> None: ...
    def resources(self) -> set[Resource]: ...
    def clear(self) -> None: ...
    def library_url(self, library: Library) -> str: ...

class DummyNeededResources:
    def need(self, resource: Resource | Group, slots: dict[Slot, Resource] | None = None) -> None: ...
    def has_resources(self) -> Literal[False]: ...

thread_local_needed_data: local

def init_needed(
    versioning: bool = False,
    versioning_use_md5: bool = False,
    recompute_hashes: bool = True,
    base_url: str | None = None,
    script_name: str | None = None,
    publisher_signature: str = ...,
    resources: Iterable[Dependable] | None = None,
) -> NeededResources: ...
def del_needed() -> None: ...
def get_needed() -> NeededResources | DummyNeededResources: ...
def clear_needed() -> None: ...

class Bundle(Renderable):
    def __init__(self) -> None: ...
    @property
    def dirname(self) -> str: ...
    @property
    def library(self) -> Library: ...
    @property
    def renderer(self) -> _Renderer: ...
    @property
    def ext(self) -> str: ...
    @property
    def relpath(self) -> str: ...
    def resources(self) -> list[Resource]: ...
    def render(self, library_url: str) -> str: ...
    def fits(self, resource: Resource) -> bool: ...
    def append(self, resource: Resource) -> None: ...
    def add_to_list(self, result: list[Renderable]) -> None: ...
