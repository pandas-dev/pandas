from _typeshed import StrPath

def augpath(
    path: StrPath,
    suffix: str = "",
    prefix: str = "",
    ext: str | None = None,
    base: str | None = None,
    dpath: str | None = None,
    multidot: bool = False,
) -> str: ...
def shrinkuser(path: StrPath, home: str = "~") -> str: ...
def expandpath(path: StrPath) -> str: ...

__all__ = ["augpath", "shrinkuser", "expandpath"]
