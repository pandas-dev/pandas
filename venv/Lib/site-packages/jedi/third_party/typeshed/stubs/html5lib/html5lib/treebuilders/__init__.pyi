from typing import Literal

treeBuilderCache: dict[str, type]

def getTreeBuilder(treeType: Literal["dom", "etree", "lxml"], implementation=None, **kwargs): ...
