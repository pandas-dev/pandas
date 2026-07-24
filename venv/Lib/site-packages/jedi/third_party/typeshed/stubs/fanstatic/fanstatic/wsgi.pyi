from _typeshed.wsgi import WSGIApplication
from typing import Any

from fanstatic.core import Resource
from fanstatic.injector import InjectorPlugin
from fanstatic.publisher import Delegator
from webob import Request, Response
from webob.dec import wsgify

def Fanstatic(
    app: WSGIApplication, publisher_signature: str = "fanstatic", injector: InjectorPlugin | None = None, **config: Any
) -> Delegator: ...
def make_fanstatic(app: WSGIApplication, global_config: Any, **local_config: Any) -> Delegator: ...

class Serf:
    resource: Resource
    def __init__(self, resource: Resource) -> None: ...
    @wsgify
    def __call__(self, request: Request) -> Response: ...

def make_serf(global_config: Any, **local_config: Any) -> Serf: ...
def resolve(name: str, module: str | None = None) -> Any: ...
