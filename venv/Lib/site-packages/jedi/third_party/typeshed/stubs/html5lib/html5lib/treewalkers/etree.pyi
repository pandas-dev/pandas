import re
from collections.abc import Callable
from types import ModuleType

tag_regexp: re.Pattern[str]

def getETreeBuilder(ElementTreeImplementation): ...

getETreeModule: Callable[..., ModuleType]
