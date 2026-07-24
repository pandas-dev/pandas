from collections import deque

from numba.core import types, cgutils



class DataPacker(object):
    """
    A helper to pack a number of typed arguments into a data structure.
    Omitted arguments (i.e. values with the type `Omitted`) are automatically
    skipped.
    """
    # XXX should DataPacker be a model for a dedicated type?

    def __init__(self, dmm, fe_types):
        self._dmm = dmm
        self._fe_types = fe_types
        self._models = [dmm.lookup(ty) for ty in fe_types]

        self._pack_map = []
        self._be_types = []
        for i, ty in enumerate(fe_types):
            if not isinstance(ty, types.Omitted):
                self._pack_map.append(i)
                self._be_types.append(self._models[i].get_data_type())

    def as_data(self, builder, values):
        """
        Return the given values packed as a data structure.
        """
        elems = [self._models[i].as_data(builder, values[i])
                 for i in self._pack_map]
        return cgutils.make_anonymous_struct(builder, elems)

    def _do_load(self, builder, ptr, formal_list=None):
        res = []
        for i, i_formal in enumerate(self._pack_map):
            elem_ptr = cgutils.gep_inbounds(builder, ptr, 0, i)
            val = self._models[i_formal].load_from_data_pointer(builder, elem_ptr)
            if formal_list is None:
                res.append((self._fe_types[i_formal], val))
            else:
                formal_list[i_formal] = val
        return res

    def load(self, builder, ptr):
        """
        Load the packed values and return a (type, value) tuples.
        """
        return self._do_load(builder, ptr)

    def load_into(self, builder, ptr, formal_list):
        """
        Load the packed values into a sequence indexed by formal
        argument number (skipping any Omitted position).
        """
        self._do_load(builder, ptr, formal_list)


class ArgPacker(object):
    """
    Compute the position for each high-level typed argument.
    It flattens every composite argument into primitive types.
    It maintains a position map for unflattening the arguments.

    Since struct (esp. nested struct) have specific ABI requirements (e.g.
    alignment, pointer address-space, ...) in different architecture (e.g.
    OpenCL, CUDA), flattening composite argument types simplifes the call
    setup from the Python side.  Functions are receiving simple primitive
    types and there are only a handful of these.
    """

    def __init__(self, dmm, fe_args):
        self._dmm = dmm
        self._fe_args = fe_args
        self._nargs = len(fe_args)

        self._dm_args = []
        argtys = []
        for ty in fe_args:
            dm = self._dmm.lookup(ty)
            self._dm_args.append(dm)
            argtys.append(dm.get_argument_type())
        self._unflattener = _Unflattener(argtys)
        self._be_args = list(_flatten(argtys))

    def as_arguments(self, builder, values):
        """Flatten all argument values
        """
        if len(values) != self._nargs:
            raise TypeError("invalid number of args: expected %d, got %d"
                            % (self._nargs, len(values)))

        if not values:
            return ()

        args = [dm.as_argument(builder, val)
                for dm, val in zip(self._dm_args, values)
                ]

        args = tuple(_flatten(args))
        return args

    def from_arguments(self, builder, args):
        """Unflatten all argument values
        """

        valtree = self._unflattener.unflatten(args)
        values = [dm.from_argument(builder, val)
                  for dm, val in zip(self._dm_args, valtree)
                  ]

        return values

    def assign_names(self, args, names):
        """Assign names for each flattened argument values.
        """

        valtree = self._unflattener.unflatten(args)
        for aval, aname in zip(valtree, names):
            self._assign_names(aval, aname)

    def _assign_names(self, val_or_nested, name, depth=()):
        if isinstance(val_or_nested, (tuple, list)):
            for pos, aval in enumerate(val_or_nested):
                self._assign_names(aval, name, depth=depth + (pos,))
        else:
            postfix = '.'.join(map(str, depth))
            parts = [name, postfix]
            val_or_nested.name = '.'.join(filter(bool, parts))

    @property
    def argument_types(self):
        """Return a list of LLVM types that are results of flattening
        composite types.
        """
        return tuple(ty for ty in self._be_args if ty != ())


def _flatten(iterable):
    """
    Flatten nested iterable of (tuple, list).
    """
    def rec(iterable):
        for i in iterable:
            if isinstance(i, (tuple, list)):
                for j in rec(i):
                    yield j
            else:
                yield i
    return rec(iterable)


_PUSH_LIST = 1
_APPEND_NEXT_VALUE = 2
_APPEND_EMPTY_TUPLE = 3
_POP = 4

class _Unflattener(object):
    """
    An object used to unflatten nested sequences after a given pattern
    (an arbitrarily nested sequence).
    The pattern shows the nested sequence shape desired when unflattening;
    the values it contains are irrelevant.
    """

    def __init__(self, pattern):
        self._code = self._build_unflatten_code(pattern)

    def _build_unflatten_code(self, iterable):
        """Build the unflatten opcode sequence for the given *iterable* structure
        (an iterable of nested sequences).
        """
        code = []
        def rec(iterable):
            for i in iterable:
                if isinstance(i, (tuple, list)):
                    if len(i) > 0:
                        code.append(_PUSH_LIST)
                        rec(i)
                        code.append(_POP)
                    else:
                        code.append(_APPEND_EMPTY_TUPLE)
                else:
                    code.append(_APPEND_NEXT_VALUE)

        rec(iterable)
        return code

    def unflatten(self, flatiter):
        """Rebuild a nested tuple structure.
        """
        vals = deque(flatiter)

        res = []
        cur = res
        stack = []
        for op in self._code:
            if op is _PUSH_LIST:
                stack.append(cur)
                cur.append([])
                cur = cur[-1]
            elif op is _APPEND_NEXT_VALUE:
                cur.append(vals.popleft())
            elif op is _APPEND_EMPTY_TUPLE:
                cur.append(())
            elif op is _POP:
                cur = stack.pop()

        assert not stack, stack
        assert not vals, vals

        return res
