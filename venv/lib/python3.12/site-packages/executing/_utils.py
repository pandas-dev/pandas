
import ast
import sys
import dis
from typing import cast, Any,Iterator
import types



def assert_(condition, message=""):
    # type: (Any, str) -> None
    """
    Like an assert statement, but unaffected by -O
    :param condition: value that is expected to be truthy
    :type message: Any
    """
    if not condition:
        raise AssertionError(str(message))


if sys.version_info >= (3, 4):
    # noinspection PyUnresolvedReferences
    _get_instructions = dis.get_instructions
    from dis import Instruction as _Instruction
    
    class Instruction(_Instruction):
        lineno = None  # type: int
else:
    from collections import namedtuple

    class Instruction(namedtuple('Instruction', 'offset argval opname starts_line')):
        lineno = None # type: int

    from dis import HAVE_ARGUMENT, EXTENDED_ARG, hasconst, opname, findlinestarts, hasname

    # Based on dis.disassemble from 2.7
    # Left as similar as possible for easy diff

    def _get_instructions(co):
        # type: (types.CodeType) -> Iterator[Instruction]
        code = co.co_code
        linestarts = dict(findlinestarts(co))
        n = len(code)
        i = 0
        extended_arg = 0
        while i < n:
            offset = i
            c = code[i]
            op = ord(c)
            lineno = linestarts.get(i)
            argval = None
            i = i + 1
            if op >= HAVE_ARGUMENT:
                oparg = ord(code[i]) + ord(code[i + 1]) * 256 + extended_arg
                extended_arg = 0
                i = i + 2
                if op == EXTENDED_ARG:
                    extended_arg = oparg * 65536

                if op in hasconst:
                    argval = co.co_consts[oparg]
                elif op in hasname:
                    argval = co.co_names[oparg]
                elif opname[op] == 'LOAD_FAST':
                    argval = co.co_varnames[oparg]
            yield Instruction(offset, argval, opname[op], lineno)

def get_instructions(co):
    # type: (types.CodeType) -> Iterator[EnhancedInstruction]
    lineno = co.co_firstlineno
    for inst in _get_instructions(co):
        inst = cast(EnhancedInstruction, inst)
        lineno = inst.starts_line or lineno
        assert_(lineno)
        inst.lineno = lineno
        yield inst


# Type class used to expand out the definition of AST to include fields added by this library
# It's not actually used for anything other than type checking though!
class EnhancedAST(ast.AST):
    parent = None  # type: EnhancedAST

# Type class used to expand out the definition of AST to include fields added by this library
# It's not actually used for anything other than type checking though!
class EnhancedInstruction(Instruction):
    _copied = None # type: bool





def mangled_name(node):
    # type: (EnhancedAST) -> str
    """

    Parameters:
        node: the node which should be mangled
        name: the name of the node

    Returns:
        The mangled name of `node`
    """

    function_class_types=(ast.FunctionDef, ast.ClassDef, ast.AsyncFunctionDef)

    if isinstance(node, ast.Attribute):
        name = node.attr
    elif isinstance(node, ast.Name):
        name = node.id
    elif isinstance(node, (ast.alias)):
        name = node.asname or node.name.split(".")[0]
    elif isinstance(node, function_class_types):
        name = node.name
    elif isinstance(node, ast.ExceptHandler):
        assert node.name
        name = node.name
    elif sys.version_info >= (3,12) and isinstance(node,ast.TypeVar):
        name=node.name
    else:
        raise TypeError("no node to mangle")

    if name.startswith("__") and not name.endswith("__"):

        parent,child=node.parent,node

        while not (isinstance(parent,ast.ClassDef) and child not in parent.bases):
            if not hasattr(parent,"parent"):
                break # pragma: no mutate

            parent,child=parent.parent,parent
        else:
            class_name=parent.name.lstrip("_")
            if class_name!="" and child not in parent.decorator_list:
                return "_" + class_name + name

            

    return name
