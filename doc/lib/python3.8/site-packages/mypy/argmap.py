"""Utilities for mapping between actual and formal arguments (and their types)."""

from typing import List, Optional, Sequence, Callable, Set

from mypy.types import (
    Type, Instance, TupleType, AnyType, TypeOfAny, TypedDictType, get_proper_type
)
from mypy import nodes


def map_actuals_to_formals(actual_kinds: List[int],
                           actual_names: Optional[Sequence[Optional[str]]],
                           formal_kinds: List[int],
                           formal_names: Sequence[Optional[str]],
                           actual_arg_type: Callable[[int],
                                                     Type]) -> List[List[int]]:
    """Calculate mapping between actual (caller) args and formals.

    The result contains a list of caller argument indexes mapping to each
    callee argument index, indexed by callee index.

    The caller_arg_type argument should evaluate to the type of the actual
    argument type with the given index.
    """
    nformals = len(formal_kinds)
    formal_to_actual = [[] for i in range(nformals)]  # type: List[List[int]]
    ambiguous_actual_kwargs = []  # type: List[int]
    fi = 0
    for ai, actual_kind in enumerate(actual_kinds):
        if actual_kind == nodes.ARG_POS:
            if fi < nformals:
                if formal_kinds[fi] in [nodes.ARG_POS, nodes.ARG_OPT,
                                        nodes.ARG_NAMED, nodes.ARG_NAMED_OPT]:
                    formal_to_actual[fi].append(ai)
                    fi += 1
                elif formal_kinds[fi] == nodes.ARG_STAR:
                    formal_to_actual[fi].append(ai)
        elif actual_kind == nodes.ARG_STAR:
            # We need to know the actual type to map varargs.
            actualt = get_proper_type(actual_arg_type(ai))
            if isinstance(actualt, TupleType):
                # A tuple actual maps to a fixed number of formals.
                for _ in range(len(actualt.items)):
                    if fi < nformals:
                        if formal_kinds[fi] != nodes.ARG_STAR2:
                            formal_to_actual[fi].append(ai)
                        else:
                            break
                        if formal_kinds[fi] != nodes.ARG_STAR:
                            fi += 1
            else:
                # Assume that it is an iterable (if it isn't, there will be
                # an error later).
                while fi < nformals:
                    if formal_kinds[fi] in (nodes.ARG_NAMED, nodes.ARG_NAMED_OPT, nodes.ARG_STAR2):
                        break
                    else:
                        formal_to_actual[fi].append(ai)
                    if formal_kinds[fi] == nodes.ARG_STAR:
                        break
                    fi += 1
        elif actual_kind in (nodes.ARG_NAMED, nodes.ARG_NAMED_OPT):
            assert actual_names is not None, "Internal error: named kinds without names given"
            name = actual_names[ai]
            if name in formal_names:
                formal_to_actual[formal_names.index(name)].append(ai)
            elif nodes.ARG_STAR2 in formal_kinds:
                formal_to_actual[formal_kinds.index(nodes.ARG_STAR2)].append(ai)
        else:
            assert actual_kind == nodes.ARG_STAR2
            actualt = get_proper_type(actual_arg_type(ai))
            if isinstance(actualt, TypedDictType):
                for name, value in actualt.items.items():
                    if name in formal_names:
                        formal_to_actual[formal_names.index(name)].append(ai)
                    elif nodes.ARG_STAR2 in formal_kinds:
                        formal_to_actual[formal_kinds.index(nodes.ARG_STAR2)].append(ai)
            else:
                # We don't exactly know which **kwargs are provided by the
                # caller, so we'll defer until all the other unambiguous
                # actuals have been processed
                ambiguous_actual_kwargs.append(ai)

    if ambiguous_actual_kwargs:
        # Assume the ambiguous kwargs will fill the remaining arguments.
        #
        # TODO: If there are also tuple varargs, we might be missing some potential
        #       matches if the tuple was short enough to not match everything.
        unmatched_formals = [fi for fi in range(nformals)
                             if (formal_names[fi]
                                 and (not formal_to_actual[fi]
                                      or actual_kinds[formal_to_actual[fi][0]] == nodes.ARG_STAR)
                                 and formal_kinds[fi] != nodes.ARG_STAR)
                             or formal_kinds[fi] == nodes.ARG_STAR2]
        for ai in ambiguous_actual_kwargs:
            for fi in unmatched_formals:
                formal_to_actual[fi].append(ai)

    return formal_to_actual


def map_formals_to_actuals(actual_kinds: List[int],
                           actual_names: Optional[Sequence[Optional[str]]],
                           formal_kinds: List[int],
                           formal_names: List[Optional[str]],
                           actual_arg_type: Callable[[int],
                                                     Type]) -> List[List[int]]:
    """Calculate the reverse mapping of map_actuals_to_formals."""
    formal_to_actual = map_actuals_to_formals(actual_kinds,
                                              actual_names,
                                              formal_kinds,
                                              formal_names,
                                              actual_arg_type)
    # Now reverse the mapping.
    actual_to_formal = [[] for _ in actual_kinds]  # type: List[List[int]]
    for formal, actuals in enumerate(formal_to_actual):
        for actual in actuals:
            actual_to_formal[actual].append(formal)
    return actual_to_formal


class ArgTypeExpander:
    """Utility class for mapping actual argument types to formal arguments.

    One of the main responsibilities is to expand caller tuple *args and TypedDict
    **kwargs, and to keep track of which tuple/TypedDict items have already been
    consumed.

    Example:

       def f(x: int, *args: str) -> None: ...
       f(*(1, 'x', 1.1))

    We'd call expand_actual_type three times:

      1. The first call would provide 'int' as the actual type of 'x' (from '1').
      2. The second call would provide 'str' as one of the actual types for '*args'.
      2. The third call would provide 'float' as one of the actual types for '*args'.

    A single instance can process all the arguments for a single call. Each call
    needs a separate instance since instances have per-call state.
    """

    def __init__(self) -> None:
        # Next tuple *args index to use.
        self.tuple_index = 0
        # Keyword arguments in TypedDict **kwargs used.
        self.kwargs_used = set()  # type: Set[str]

    def expand_actual_type(self,
                           actual_type: Type,
                           actual_kind: int,
                           formal_name: Optional[str],
                           formal_kind: int) -> Type:
        """Return the actual (caller) type(s) of a formal argument with the given kinds.

        If the actual argument is a tuple *args, return the next individual tuple item that
        maps to the formal arg.

        If the actual argument is a TypedDict **kwargs, return the next matching typed dict
        value type based on formal argument name and kind.

        This is supposed to be called for each formal, in order. Call multiple times per
        formal if multiple actuals map to a formal.
        """
        actual_type = get_proper_type(actual_type)
        if actual_kind == nodes.ARG_STAR:
            if isinstance(actual_type, Instance):
                if actual_type.type.fullname == 'builtins.list':
                    # List *arg.
                    return actual_type.args[0]
                elif actual_type.args:
                    # TODO: Try to map type arguments to Iterable
                    return actual_type.args[0]
                else:
                    return AnyType(TypeOfAny.from_error)
            elif isinstance(actual_type, TupleType):
                # Get the next tuple item of a tuple *arg.
                if self.tuple_index >= len(actual_type.items):
                    # Exhausted a tuple -- continue to the next *args.
                    self.tuple_index = 1
                else:
                    self.tuple_index += 1
                return actual_type.items[self.tuple_index - 1]
            else:
                return AnyType(TypeOfAny.from_error)
        elif actual_kind == nodes.ARG_STAR2:
            if isinstance(actual_type, TypedDictType):
                if formal_kind != nodes.ARG_STAR2 and formal_name in actual_type.items:
                    # Lookup type based on keyword argument name.
                    assert formal_name is not None
                else:
                    # Pick an arbitrary item if no specified keyword is expected.
                    formal_name = (set(actual_type.items.keys()) - self.kwargs_used).pop()
                self.kwargs_used.add(formal_name)
                return actual_type.items[formal_name]
            elif (isinstance(actual_type, Instance)
                  and (actual_type.type.fullname == 'builtins.dict')):
                # Dict **arg.
                # TODO: Handle arbitrary Mapping
                return actual_type.args[1]
            else:
                return AnyType(TypeOfAny.from_error)
        else:
            # No translation for other kinds -- 1:1 mapping.
            return actual_type
