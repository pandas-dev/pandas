from _typeshed.wsgi import StartResponse, WSGIApplication, WSGIEnvironment
from abc import abstractmethod
from collections.abc import Iterable
from typing import Any, Literal, TypedDict, type_check_only
from typing_extensions import Unpack

from fanstatic.core import Dependable, NeededResources, Resource
from fanstatic.inclusion import Inclusion
from webob import Request, Response

@type_check_only
class _NeededResourcesConfig(TypedDict, total=False):
    versioning: bool
    versioning_use_md5: bool
    recompute_hashes: bool
    base_url: str | None
    script_name: str | None
    publisher_signature: str
    resources: Iterable[Dependable] | None

@type_check_only
class _InjectorPluginOptions(TypedDict, total=False):
    compile: bool
    bundle: bool
    rollup: bool
    debug: bool
    minified: bool

@type_check_only
class _TopBottomInjectorPluginOptions(_InjectorPluginOptions, total=False):
    bottom: bool
    force_bottom: bool

CONTENT_TYPES: list[str]

class Injector:
    app: WSGIApplication
    config: _NeededResourcesConfig
    injector: InjectorPlugin
    def __init__(
        self, app: WSGIApplication, injector: InjectorPlugin | None = None, **config: Unpack[_NeededResourcesConfig]
    ) -> None: ...
    def __call__(self, environ: WSGIEnvironment, start_response: StartResponse) -> Iterable[bytes]: ...

class InjectorPlugin:
    @property
    @abstractmethod
    def name(self) -> str: ...
    def __init__(self, options: _InjectorPluginOptions) -> None: ...
    def make_inclusion(self, needed: NeededResources, resources: set[Resource] | None = None) -> Inclusion: ...
    def __call__(
        self, html: bytes, needed: NeededResources, request: Request | None = None, response: Response | None = None
    ) -> None: ...

class TopBottomInjector(InjectorPlugin):
    name: Literal["topbottom"]
    def __init__(self, options: _TopBottomInjectorPluginOptions) -> None: ...
    def group(self, needed: NeededResources) -> tuple[Inclusion, Inclusion]: ...

def make_injector(app: WSGIApplication, global_config: Any, **local_config: Any) -> Injector: ...
