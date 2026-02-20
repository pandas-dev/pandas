"""
    Module of utilities for transforming a lark.Tree into a custom Abstract Syntax Tree (AST defined in classes)
"""

import inspect, re
import types
from typing import Optional, Callable

from lark import Transformer, v_args

class Ast:
    """Abstract class

    Subclasses will be collected by `create_transformer()`
    """
    pass

class AsList:
    """Abstract class

    Subclasses will be instantiated with the parse results as a single list, instead of as arguments.
    """

class WithMeta:
    """Abstract class

    Subclasses will be instantiated with the Meta instance of the tree. (see ``v_args`` for more detail)
    """
    pass

def camel_to_snake(name):
    return re.sub(r'(?<!^)(?=[A-Z])', '_', name).lower()

def create_transformer(ast_module: types.ModuleType,
                       transformer: Optional[Transformer]=None,
                       decorator_factory: Callable=v_args) -> Transformer:
    """Collects `Ast` subclasses from the given module, and creates a Lark transformer that builds the AST.

    For each class, we create a corresponding rule in the transformer, with a matching name.
    CamelCase names will be converted into snake_case. Example: "CodeBlock" -> "code_block".

    Classes starting with an underscore (`_`) will be skipped.

    Parameters:
        ast_module: A Python module containing all the subclasses of ``ast_utils.Ast``
        transformer (Optional[Transformer]): An initial transformer. Its attributes may be overwritten.
        decorator_factory (Callable): An optional callable accepting two booleans, inline, and meta,
            and returning a decorator for the methods of ``transformer``. (default: ``v_args``).
    """
    t = transformer or Transformer()

    for name, obj in inspect.getmembers(ast_module):
        if not name.startswith('_') and inspect.isclass(obj):
            if issubclass(obj, Ast):
                wrapper = decorator_factory(inline=not issubclass(obj, AsList), meta=issubclass(obj, WithMeta))
                obj = wrapper(obj).__get__(t)
                setattr(t, camel_to_snake(name), obj)

    return t
