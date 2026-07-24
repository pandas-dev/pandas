import re
from _typeshed import Incomplete
from collections.abc import Callable
from types import ModuleType

tag_regexp: re.Pattern[str]

def getETreeBuilder(ElementTreeImplementation, fullTree: bool = False) -> dict[str, Incomplete]: ...

getETreeModule: Callable[..., ModuleType]
