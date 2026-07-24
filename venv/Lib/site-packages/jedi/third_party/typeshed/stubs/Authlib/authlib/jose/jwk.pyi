from _typeshed import Incomplete
from typing_extensions import deprecated

@deprecated("Please use `JsonWebKey` directly.")
def loads(obj, kid=None): ...
@deprecated("Please use `JsonWebKey` directly.")
def dumps(key, kty=None, **params) -> dict[Incomplete, Incomplete]: ...
