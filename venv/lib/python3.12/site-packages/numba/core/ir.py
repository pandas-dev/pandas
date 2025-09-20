from collections import defaultdict
import copy
import itertools
import os
import linecache
import pprint
import re
import sys
import operator
from types import FunctionType, BuiltinFunctionType
from functools import total_ordering
from io import StringIO

from numba.core import errors, config
from numba.core.utils import (BINOPS_TO_OPERATORS, INPLACE_BINOPS_TO_OPERATORS,
                              UNARY_BUITINS_TO_OPERATORS, OPERATORS_TO_BUILTINS)
from numba.core.errors import (NotDefinedError, RedefinedError,
                               VerificationError, ConstantInferenceError)
from numba.core import consts

# terminal color markup
_termcolor = errors.termcolor()


class Loc(object):
    """Source location

    """
    _defmatcher = re.compile(r'def\s+(\w+)')

    def __init__(self, filename, line, col=None, maybe_decorator=False):
        """ Arguments:
        filename - name of the file
        line - line in file
        col - column
        maybe_decorator - Set to True if location is likely a jit decorator
        """
        self.filename = filename
        self.line = line
        self.col = col
        self.lines = None # the source lines from the linecache
        self.maybe_decorator = maybe_decorator

    def __eq__(self, other):
        # equivalence is solely based on filename, line and col
        if type(self) is not type(other): return False
        if self.filename != other.filename: return False
        if self.line != other.line: return False
        if self.col != other.col: return False
        return True

    def __ne__(self, other):
        return not self.__eq__(other)

    @classmethod
    def from_function_id(cls, func_id):
        return cls(func_id.filename, func_id.firstlineno, maybe_decorator=True)

    def __repr__(self):
        return "Loc(filename=%s, line=%s, col=%s)" % (self.filename,
                                                      self.line, self.col)

    def __str__(self):
        if self.col is not None:
            return "%s (%s:%s)" % (self.filename, self.line, self.col)
        else:
            return "%s (%s)" % (self.filename, self.line)

    def _find_definition(self):
        # try and find a def, go backwards from error line
        fn_name = None
        lines = self.get_lines()
        for x in reversed(lines[:self.line - 1]):
            # the strip and startswith is to handle user code with commented out
            # 'def' or use of 'def' in a docstring.
            if x.strip().startswith('def '):
                fn_name = x
                break

        return fn_name

    def _raw_function_name(self):
        defn = self._find_definition()
        if defn:
            m = self._defmatcher.match(defn.strip())
            if m:
                return m.groups()[0]
        # Probably exec(<string>) or REPL.
        return None

    def get_lines(self):
        if self.lines is None:
            path = self._get_path()
            # Avoid reading from dynamic string. They are most likely
            # overridden. Problem started with Python 3.13. "<string>" seems
            # to be something from multiprocessing.
            lns = [] if path == "<string>" else linecache.getlines(path)
            self.lines = lns
        return self.lines

    def _get_path(self):
        path = None
        try:
            # Try to get a relative path
            # ipython/jupyter input just returns as self.filename
            path = os.path.relpath(self.filename)
        except ValueError:
            # Fallback to absolute path if error occurred in getting the
            # relative path.
            # This may happen on windows if the drive is different
            path = os.path.abspath(self.filename)
        return path


    def strformat(self, nlines_up=2):

        lines = self.get_lines()

        use_line = self.line

        if self.maybe_decorator:
            # try and sort out a better `loc`, if it's suspected that this loc
            # points at a jit decorator by virtue of
            # `__code__.co_firstlineno`

            # get lines, add a dummy entry at the start as lines count from
            # 1 but list index counts from 0
            tmplines = [''] + lines

            if lines and use_line and 'def ' not in tmplines[use_line]:
                # look forward 10 lines, unlikely anyone managed to stretch
                # a jit call declaration over >10 lines?!
                min_line = max(0, use_line)
                max_line = use_line + 10
                selected = tmplines[min_line : max_line]
                index = 0
                for idx, x in enumerate(selected):
                    if 'def ' in x:
                        index = idx
                        break
                use_line = use_line + index


        ret = [] # accumulates output
        if lines and use_line > 0:

            def count_spaces(string):
                spaces = 0
                for x in itertools.takewhile(str.isspace, str(string)):
                    spaces += 1
                return spaces

            # A few places in the code still use no `loc` or default to line 1
            # this is often in places where exceptions are used for the purposes
            # of flow control. As a result max is in use to prevent slice from
            # `[negative: positive]`
            selected = lines[max(0, use_line - nlines_up):use_line]

            # see if selected contains a definition
            def_found = False
            for x in selected:
                if 'def ' in x:
                    def_found = True

            # no definition found, try and find one
            if not def_found:
                # try and find a def, go backwards from error line
                fn_name = None
                for x in reversed(lines[:use_line - 1]):
                    if 'def ' in x:
                        fn_name = x
                        break
                if fn_name:
                    ret.append(fn_name)
                    spaces = count_spaces(x)
                    ret.append(' '*(4 + spaces) + '<source elided>\n')

            if selected:
                ret.extend(selected[:-1])
                ret.append(_termcolor.highlight(selected[-1]))

                # point at the problem with a caret
                spaces = count_spaces(selected[-1])
                ret.append(' '*(spaces) + _termcolor.indicate("^"))

        # if in the REPL source may not be available
        if not ret:
            if not lines:
                ret = "<source missing, REPL/exec in use?>"
            elif use_line <= 0:
                ret = "<source line number missing>"


        err = _termcolor.filename('\nFile "%s", line %d:')+'\n%s'
        tmp = err % (self._get_path(), use_line, _termcolor.code(''.join(ret)))
        return tmp

    def with_lineno(self, line, col=None):
        """
        Return a new Loc with this line number.
        """
        return type(self)(self.filename, line, col)

    def short(self):
        """
        Returns a short string
        """
        shortfilename = os.path.basename(self.filename)
        return "%s:%s" % (shortfilename, self.line)


# Used for annotating errors when source location is unknown.
unknown_loc = Loc("unknown location", 0, 0)


@total_ordering
class SlotEqualityCheckMixin(object):
    # some ir nodes are __dict__ free using __slots__ instead, this mixin
    # should not trigger the unintended creation of __dict__.
    __slots__ = tuple()

    def __eq__(self, other):
        if type(self) is type(other):
            for name in self.__slots__:
                if getattr(self, name) != getattr(other, name):
                    return False
            else:
                return True
        return False

    def __le__(self, other):
        return str(self) <= str(other)

    def __hash__(self):
        return id(self)


@total_ordering
class EqualityCheckMixin(object):
    """ Mixin for basic equality checking """

    def __eq__(self, other):
        if type(self) is type(other):
            def fixup(adict):
                bad = ('loc', 'scope')
                d = dict(adict)
                for x in bad:
                    d.pop(x, None)
                return d
            d1 = fixup(self.__dict__)
            d2 = fixup(other.__dict__)
            if d1 == d2:
                return True
        return False

    def __le__(self, other):
        return str(self) < str(other)

    def __hash__(self):
        return id(self)


class VarMap(object):
    def __init__(self):
        self._con = {}

    def define(self, name, var):
        if name in self._con:
            raise RedefinedError(name)
        else:
            self._con[name] = var

    def get(self, name):
        try:
            return self._con[name]
        except KeyError:
            raise NotDefinedError(name)

    def __contains__(self, name):
        return name in self._con

    def __len__(self):
        return len(self._con)

    def __repr__(self):
        return pprint.pformat(self._con)

    def __hash__(self):
        return hash(self.name)

    def __iter__(self):
        return self._con.iterkeys()

    def __eq__(self, other):
        if type(self) is type(other):
            # check keys only, else __eq__ ref cycles, scope -> varmap -> var
            return self._con.keys() == other._con.keys()
        return False

    def __ne__(self, other):
        return not self.__eq__(other)


class AbstractRHS(object):
    """Abstract base class for anything that can be the RHS of an assignment.
    This class **does not** define any methods.
    """


class Inst(EqualityCheckMixin, AbstractRHS):
    """
    Base class for all IR instructions.
    """

    def list_vars(self):
        """
        List the variables used (read or written) by the instruction.
        """
        raise NotImplementedError

    def _rec_list_vars(self, val):
        """
        A recursive helper used to implement list_vars() in subclasses.
        """
        if isinstance(val, Var):
            return [val]
        elif isinstance(val, Inst):
            return val.list_vars()
        elif isinstance(val, (list, tuple)):
            lst = []
            for v in val:
                lst.extend(self._rec_list_vars(v))
            return lst
        elif isinstance(val, dict):
            lst = []
            for v in val.values():
                lst.extend(self._rec_list_vars(v))
            return lst
        else:
            return []


class Stmt(Inst):
    """
    Base class for IR statements (instructions which can appear on their
    own in a Block).
    """
    # Whether this statement ends its basic block (i.e. it will either jump
    # to another block or exit the function).
    is_terminator = False
    # Whether this statement exits the function.
    is_exit = False

    def list_vars(self):
        return self._rec_list_vars(self.__dict__)


class Terminator(Stmt):
    """
    IR statements that are terminators: the last statement in a block.
    A terminator must either:
    - exit the function
    - jump to a block

    All subclass of Terminator must override `.get_targets()` to return a list
    of jump targets.
    """
    is_terminator = True

    def get_targets(self):
        raise NotImplementedError(type(self))


class Expr(Inst):
    """
    An IR expression (an instruction which can only be part of a larger
    statement).
    """

    def __init__(self, op, loc, **kws):
        assert isinstance(op, str)
        assert isinstance(loc, Loc)
        self.op = op
        self.loc = loc
        self._kws = kws

    def __getattr__(self, name):
        if name.startswith('_'):
            return Inst.__getattr__(self, name)
        return self._kws[name]

    def __setattr__(self, name, value):
        if name in ('op', 'loc', '_kws'):
            self.__dict__[name] = value
        else:
            self._kws[name] = value

    @classmethod
    def binop(cls, fn, lhs, rhs, loc):
        assert isinstance(fn, BuiltinFunctionType)
        assert isinstance(lhs, Var)
        assert isinstance(rhs, Var)
        assert isinstance(loc, Loc)
        op = 'binop'
        return cls(op=op, loc=loc, fn=fn, lhs=lhs, rhs=rhs,
                   static_lhs=UNDEFINED, static_rhs=UNDEFINED)

    @classmethod
    def inplace_binop(cls, fn, immutable_fn, lhs, rhs, loc):
        assert isinstance(fn, BuiltinFunctionType)
        assert isinstance(immutable_fn, BuiltinFunctionType)
        assert isinstance(lhs, Var)
        assert isinstance(rhs, Var)
        assert isinstance(loc, Loc)
        op = 'inplace_binop'
        return cls(op=op, loc=loc, fn=fn, immutable_fn=immutable_fn,
                   lhs=lhs, rhs=rhs,
                   static_lhs=UNDEFINED, static_rhs=UNDEFINED)

    @classmethod
    def unary(cls, fn, value, loc):
        assert isinstance(value, (str, Var, FunctionType))
        assert isinstance(loc, Loc)
        op = 'unary'
        fn = UNARY_BUITINS_TO_OPERATORS.get(fn, fn)
        return cls(op=op, loc=loc, fn=fn, value=value)

    @classmethod
    def call(cls, func, args, kws, loc, vararg=None, varkwarg=None, target=None):
        assert isinstance(func, Var)
        assert isinstance(loc, Loc)
        op = 'call'
        return cls(op=op, loc=loc, func=func, args=args, kws=kws,
                   vararg=vararg, varkwarg=varkwarg, target=target)

    @classmethod
    def build_tuple(cls, items, loc):
        assert isinstance(loc, Loc)
        op = 'build_tuple'
        return cls(op=op, loc=loc, items=items)

    @classmethod
    def build_list(cls, items, loc):
        assert isinstance(loc, Loc)
        op = 'build_list'
        return cls(op=op, loc=loc, items=items)

    @classmethod
    def build_set(cls, items, loc):
        assert isinstance(loc, Loc)
        op = 'build_set'
        return cls(op=op, loc=loc, items=items)

    @classmethod
    def build_map(cls, items, size, literal_value, value_indexes, loc):
        assert isinstance(loc, Loc)
        op = 'build_map'
        return cls(op=op, loc=loc, items=items, size=size,
                   literal_value=literal_value, value_indexes=value_indexes)

    @classmethod
    def pair_first(cls, value, loc):
        assert isinstance(value, Var)
        op = 'pair_first'
        return cls(op=op, loc=loc, value=value)

    @classmethod
    def pair_second(cls, value, loc):
        assert isinstance(value, Var)
        assert isinstance(loc, Loc)
        op = 'pair_second'
        return cls(op=op, loc=loc, value=value)

    @classmethod
    def getiter(cls, value, loc):
        assert isinstance(value, Var)
        assert isinstance(loc, Loc)
        op = 'getiter'
        return cls(op=op, loc=loc, value=value)

    @classmethod
    def iternext(cls, value, loc):
        assert isinstance(value, Var)
        assert isinstance(loc, Loc)
        op = 'iternext'
        return cls(op=op, loc=loc, value=value)

    @classmethod
    def exhaust_iter(cls, value, count, loc):
        assert isinstance(value, Var)
        assert isinstance(count, int)
        assert isinstance(loc, Loc)
        op = 'exhaust_iter'
        return cls(op=op, loc=loc, value=value, count=count)

    @classmethod
    def getattr(cls, value, attr, loc):
        assert isinstance(value, Var)
        assert isinstance(attr, str)
        assert isinstance(loc, Loc)
        op = 'getattr'
        return cls(op=op, loc=loc, value=value, attr=attr)

    @classmethod
    def getitem(cls, value, index, loc):
        assert isinstance(value, Var)
        assert isinstance(index, Var)
        assert isinstance(loc, Loc)
        op = 'getitem'
        fn = operator.getitem
        return cls(op=op, loc=loc, value=value, index=index, fn=fn)

    @classmethod
    def typed_getitem(cls, value, dtype, index, loc):
        assert isinstance(value, Var)
        assert isinstance(loc, Loc)
        op = 'typed_getitem'
        return cls(op=op, loc=loc, value=value, dtype=dtype,
                   index=index)

    @classmethod
    def static_getitem(cls, value, index, index_var, loc):
        assert isinstance(value, Var)
        assert index_var is None or isinstance(index_var, Var)
        assert isinstance(loc, Loc)
        op = 'static_getitem'
        fn = operator.getitem
        return cls(op=op, loc=loc, value=value, index=index,
                   index_var=index_var, fn=fn)

    @classmethod
    def cast(cls, value, loc):
        """
        A node for implicit casting at the return statement
        """
        assert isinstance(value, Var)
        assert isinstance(loc, Loc)
        op = 'cast'
        return cls(op=op, value=value, loc=loc)

    @classmethod
    def phi(cls, loc):
        """Phi node
        """
        assert isinstance(loc, Loc)
        return cls(op='phi', incoming_values=[], incoming_blocks=[], loc=loc)

    @classmethod
    def make_function(cls, name, code, closure, defaults, loc):
        """
        A node for making a function object.
        """
        assert isinstance(loc, Loc)
        op = 'make_function'
        return cls(op=op, name=name, code=code, closure=closure, defaults=defaults, loc=loc)

    @classmethod
    def null(cls, loc):
        """
        A node for null value.

        This node is not handled by type inference. It is only added by
        post-typing passes.
        """
        assert isinstance(loc, Loc)
        op = 'null'
        return cls(op=op, loc=loc)

    @classmethod
    def undef(cls, loc):
        """
        A node for undefined value specifically from LOAD_FAST_AND_CLEAR opcode.
        """
        assert isinstance(loc, Loc)
        op = 'undef'
        return cls(op=op, loc=loc)

    @classmethod
    def dummy(cls, op, info, loc):
        """
        A node for a dummy value.

        This node is a place holder for carrying information through to a point
        where it is rewritten into something valid. This node is not handled
        by type inference or lowering. It's presence outside of the interpreter
        renders IR as illegal.
        """
        assert isinstance(loc, Loc)
        assert isinstance(op, str)
        return cls(op=op, info=info, loc=loc)

    def __repr__(self):
        if self.op == 'call':
            args = ', '.join(str(a) for a in self.args)
            pres_order = self._kws.items() if config.DIFF_IR == 0 else sorted(self._kws.items())
            kws = ', '.join('%s=%s' % (k, v) for k, v in pres_order)
            vararg = '*%s' % (self.vararg,) if self.vararg is not None else ''
            arglist = ', '.join(filter(None, [args, vararg, kws]))
            return 'call %s(%s)' % (self.func, arglist)
        elif self.op == 'binop':
            lhs, rhs = self.lhs, self.rhs
            if self.fn == operator.contains:
                lhs, rhs = rhs, lhs
            fn = OPERATORS_TO_BUILTINS.get(self.fn, self.fn)
            return '%s %s %s' % (lhs, fn, rhs)
        else:
            pres_order = self._kws.items() if config.DIFF_IR == 0 else sorted(self._kws.items())
            args = ('%s=%s' % (k, v) for k, v in pres_order)
            return '%s(%s)' % (self.op, ', '.join(args))

    def list_vars(self):
        return self._rec_list_vars(self._kws)

    def infer_constant(self):
        raise ConstantInferenceError('%s' % self, loc=self.loc)


class SetItem(Stmt):
    """
    target[index] = value
    """

    def __init__(self, target, index, value, loc):
        assert isinstance(target, Var)
        assert isinstance(index, Var)
        assert isinstance(value, Var)
        assert isinstance(loc, Loc)
        self.target = target
        self.index = index
        self.value = value
        self.loc = loc

    def __repr__(self):
        return '%s[%s] = %s' % (self.target, self.index, self.value)


class StaticSetItem(Stmt):
    """
    target[constant index] = value
    """

    def __init__(self, target, index, index_var, value, loc):
        assert isinstance(target, Var)
        assert not isinstance(index, Var)
        assert isinstance(index_var, Var)
        assert isinstance(value, Var)
        assert isinstance(loc, Loc)
        self.target = target
        self.index = index
        self.index_var = index_var
        self.value = value
        self.loc = loc

    def __repr__(self):
        return '%s[%r] = %s' % (self.target, self.index, self.value)


class DelItem(Stmt):
    """
    del target[index]
    """

    def __init__(self, target, index, loc):
        assert isinstance(target, Var)
        assert isinstance(index, Var)
        assert isinstance(loc, Loc)
        self.target = target
        self.index = index
        self.loc = loc

    def __repr__(self):
        return 'del %s[%s]' % (self.target, self.index)


class SetAttr(Stmt):
    def __init__(self, target, attr, value, loc):
        assert isinstance(target, Var)
        assert isinstance(attr, str)
        assert isinstance(value, Var)
        assert isinstance(loc, Loc)
        self.target = target
        self.attr = attr
        self.value = value
        self.loc = loc

    def __repr__(self):
        return '(%s).%s = %s' % (self.target, self.attr, self.value)


class DelAttr(Stmt):
    def __init__(self, target, attr, loc):
        assert isinstance(target, Var)
        assert isinstance(attr, str)
        assert isinstance(loc, Loc)
        self.target = target
        self.attr = attr
        self.loc = loc

    def __repr__(self):
        return 'del (%s).%s' % (self.target, self.attr)


class StoreMap(Stmt):
    def __init__(self, dct, key, value, loc):
        assert isinstance(dct, Var)
        assert isinstance(key, Var)
        assert isinstance(value, Var)
        assert isinstance(loc, Loc)
        self.dct = dct
        self.key = key
        self.value = value
        self.loc = loc

    def __repr__(self):
        return '%s[%s] = %s' % (self.dct, self.key, self.value)


class Del(Stmt):
    def __init__(self, value, loc):
        assert isinstance(value, str)
        assert isinstance(loc, Loc)
        self.value = value
        self.loc = loc

    def __str__(self):
        return "del %s" % self.value


class Raise(Terminator):
    is_exit = True

    def __init__(self, exception, loc):
        assert exception is None or isinstance(exception, Var)
        assert isinstance(loc, Loc)
        self.exception = exception
        self.loc = loc

    def __str__(self):
        return "raise %s" % self.exception

    def get_targets(self):
        return []


class StaticRaise(Terminator):
    """
    Raise an exception class and arguments known at compile-time.
    Note that if *exc_class* is None, a bare "raise" statement is implied
    (i.e. re-raise the current exception).
    """
    is_exit = True

    def __init__(self, exc_class, exc_args, loc):
        assert exc_class is None or isinstance(exc_class, type)
        assert isinstance(loc, Loc)
        assert exc_args is None or isinstance(exc_args, tuple)
        self.exc_class = exc_class
        self.exc_args = exc_args
        self.loc = loc

    def __str__(self):
        if self.exc_class is None:
            return "<static> raise"
        elif self.exc_args is None:
            return "<static> raise %s" % (self.exc_class,)
        else:
            return "<static> raise %s(%s)" % (self.exc_class,
                                     ", ".join(map(repr, self.exc_args)))

    def get_targets(self):
        return []


class DynamicRaise(Terminator):
    """
    Raise an exception class and some argument *values* unknown at compile-time.
    Note that if *exc_class* is None, a bare "raise" statement is implied
    (i.e. re-raise the current exception).
    """
    is_exit = True

    def __init__(self, exc_class, exc_args, loc):
        assert exc_class is None or isinstance(exc_class, type)
        assert isinstance(loc, Loc)
        assert exc_args is None or isinstance(exc_args, tuple)
        self.exc_class = exc_class
        self.exc_args = exc_args
        self.loc = loc

    def __str__(self):
        if self.exc_class is None:
            return "<dynamic> raise"
        elif self.exc_args is None:
            return "<dynamic> raise %s" % (self.exc_class,)
        else:
            return "<dynamic> raise %s(%s)" % (self.exc_class,
                                     ", ".join(map(repr, self.exc_args)))

    def get_targets(self):
        return []


class TryRaise(Stmt):
    """A raise statement inside a try-block
    Similar to ``Raise`` but does not terminate.
    """
    def __init__(self, exception, loc):
        assert exception is None or isinstance(exception, Var)
        assert isinstance(loc, Loc)
        self.exception = exception
        self.loc = loc

    def __str__(self):
        return "try_raise %s" % self.exception


class StaticTryRaise(Stmt):
    """A raise statement inside a try-block.
    Similar to ``StaticRaise`` but does not terminate.
    """
    def __init__(self, exc_class, exc_args, loc):
        assert exc_class is None or isinstance(exc_class, type)
        assert isinstance(loc, Loc)
        assert exc_args is None or isinstance(exc_args, tuple)
        self.exc_class = exc_class
        self.exc_args = exc_args
        self.loc = loc

    def __str__(self):
        if self.exc_class is None:
            return f"static_try_raise"
        elif self.exc_args is None:
            return f"static_try_raise {self.exc_class}"
        else:
            args = ", ".join(map(repr, self.exc_args))
            return f"static_try_raise {self.exc_class}({args})"


class DynamicTryRaise(Stmt):
    """A raise statement inside a try-block.
    Similar to ``DynamicRaise`` but does not terminate.
    """
    def __init__(self, exc_class, exc_args, loc):
        assert exc_class is None or isinstance(exc_class, type)
        assert isinstance(loc, Loc)
        assert exc_args is None or isinstance(exc_args, tuple)
        self.exc_class = exc_class
        self.exc_args = exc_args
        self.loc = loc

    def __str__(self):
        if self.exc_class is None:
            return f"dynamic_try_raise"
        elif self.exc_args is None:
            return f"dynamic_try_raise {self.exc_class}"
        else:
            args = ", ".join(map(repr, self.exc_args))
            return f"dynamic_try_raise {self.exc_class}({args})"


class Return(Terminator):
    """
    Return to caller.
    """
    is_exit = True

    def __init__(self, value, loc):
        assert isinstance(value, Var), type(value)
        assert isinstance(loc, Loc)
        self.value = value
        self.loc = loc

    def __str__(self):
        return 'return %s' % self.value

    def get_targets(self):
        return []


class Jump(Terminator):
    """
    Unconditional branch.
    """

    def __init__(self, target, loc):
        assert isinstance(loc, Loc)
        self.target = target
        self.loc = loc

    def __str__(self):
        return 'jump %s' % self.target

    def get_targets(self):
        return [self.target]


class Branch(Terminator):
    """
    Conditional branch.
    """

    def __init__(self, cond, truebr, falsebr, loc):
        assert isinstance(cond, Var)
        assert isinstance(loc, Loc)
        self.cond = cond
        self.truebr = truebr
        self.falsebr = falsebr
        self.loc = loc

    def __str__(self):
        return 'branch %s, %s, %s' % (self.cond, self.truebr, self.falsebr)

    def get_targets(self):
        return [self.truebr, self.falsebr]


class Assign(Stmt):
    """
    Assign to a variable.
    """
    def __init__(self, value, target, loc):
        assert isinstance(value, AbstractRHS)
        assert isinstance(target, Var)
        assert isinstance(loc, Loc)
        self.value = value
        self.target = target
        self.loc = loc

    def __str__(self):
        return '%s = %s' % (self.target, self.value)


class Print(Stmt):
    """
    Print some values.
    """
    def __init__(self, args, vararg, loc):
        assert all(isinstance(x, Var) for x in args)
        assert vararg is None or isinstance(vararg, Var)
        assert isinstance(loc, Loc)
        self.args = tuple(args)
        self.vararg = vararg
        # Constant-inferred arguments
        self.consts = {}
        self.loc = loc

    def __str__(self):
        return 'print(%s)' % ', '.join(str(v) for v in self.args)


class Yield(Inst):
    def __init__(self, value, loc, index):
        assert isinstance(value, Var)
        assert isinstance(loc, Loc)
        self.value = value
        self.loc = loc
        self.index = index

    def __str__(self):
        return 'yield %s' % (self.value,)

    def list_vars(self):
        return [self.value]


class EnterWith(Stmt):
    """Enter a "with" context
    """
    def __init__(self, contextmanager, begin, end, loc):
        """
        Parameters
        ----------
        contextmanager : IR value
        begin, end : int
            The beginning and the ending offset of the with-body.
        loc : ir.Loc instance
            Source location
        """
        assert isinstance(contextmanager, Var)
        assert isinstance(loc, Loc)
        self.contextmanager = contextmanager
        self.begin = begin
        self.end = end
        self.loc = loc

    def __str__(self):
        return 'enter_with {}'.format(self.contextmanager)

    def list_vars(self):
        return [self.contextmanager]


class PopBlock(Stmt):
    """Marker statement for a pop block op code"""
    def __init__(self, loc):
        assert isinstance(loc, Loc)
        self.loc = loc

    def __str__(self):
        return 'pop_block'


class Arg(EqualityCheckMixin, AbstractRHS):
    def __init__(self, name, index, loc):
        assert isinstance(name, str)
        assert isinstance(index, int)
        assert isinstance(loc, Loc)
        self.name = name
        self.index = index
        self.loc = loc

    def __repr__(self):
        return 'arg(%d, name=%s)' % (self.index, self.name)

    def infer_constant(self):
        raise ConstantInferenceError('%s' % self, loc=self.loc)


class Const(EqualityCheckMixin, AbstractRHS):
    def __init__(self, value, loc, use_literal_type=True):
        assert isinstance(loc, Loc)
        self.value = value
        self.loc = loc
        # Note: need better way to tell if this is a literal or not.
        self.use_literal_type = use_literal_type

    def __repr__(self):
        return 'const(%s, %s)' % (type(self.value).__name__, self.value)

    def infer_constant(self):
        return self.value

    def __deepcopy__(self, memo):
        # Override to not copy constant values in code
        return Const(
            value=self.value, loc=self.loc,
            use_literal_type=self.use_literal_type,
        )


class Global(EqualityCheckMixin, AbstractRHS):
    def __init__(self, name, value, loc):
        assert isinstance(loc, Loc)
        self.name = name
        self.value = value
        self.loc = loc

    def __str__(self):
        return 'global(%s: %s)' % (self.name, self.value)

    def infer_constant(self):
        return self.value

    def __deepcopy__(self, memo):
        # don't copy value since it can fail (e.g. modules)
        # value is readonly and doesn't need copying
        return Global(self.name, self.value, copy.deepcopy(self.loc))


class FreeVar(EqualityCheckMixin, AbstractRHS):
    """
    A freevar, as loaded by LOAD_DECREF.
    (i.e. a variable defined in an enclosing non-global scope)
    """

    def __init__(self, index, name, value, loc):
        assert isinstance(index, int)
        assert isinstance(name, str)
        assert isinstance(loc, Loc)
        # index inside __code__.co_freevars
        self.index = index
        # variable name
        self.name = name
        # frozen value
        self.value = value
        self.loc = loc

    def __str__(self):
        return 'freevar(%s: %s)' % (self.name, self.value)

    def infer_constant(self):
        return self.value

    def __deepcopy__(self, memo):
        # Override to not copy constant values in code
        return FreeVar(index=self.index, name=self.name, value=self.value,
                       loc=self.loc)



class Var(EqualityCheckMixin, AbstractRHS):
    """
    Attributes
    -----------
    - scope: Scope

    - name: str

    - loc: Loc
        Definition location
    """

    def __init__(self, scope, name, loc):
        # NOTE: Use of scope=None should be removed.
        assert scope is None or isinstance(scope, Scope)
        assert isinstance(name, str)
        assert isinstance(loc, Loc)
        self.scope = scope
        self.name = name
        self.loc = loc

    def __repr__(self):
        return 'Var(%s, %s)' % (self.name, self.loc.short())

    def __str__(self):
        return self.name

    @property
    def is_temp(self):
        return self.name.startswith("$")

    @property
    def unversioned_name(self):
        """The unversioned name of this variable, i.e. SSA renaming removed
        """
        for k, redef_set in self.scope.var_redefinitions.items():
            if self.name in redef_set:
                return k
        return self.name

    @property
    def versioned_names(self):
        """Known versioned names for this variable, i.e. known variable names in
        the scope that have been formed from applying SSA to this variable
        """
        return self.scope.get_versions_of(self.unversioned_name)

    @property
    def all_names(self):
        """All known versioned and unversioned names for this variable
        """
        return self.versioned_names | {self.unversioned_name,}

    def __deepcopy__(self, memo):
        out = Var(copy.deepcopy(self.scope, memo), self.name, self.loc)
        memo[id(self)] = out
        return out


class Scope(EqualityCheckMixin):
    """
    Attributes
    -----------
    - parent: Scope
        Parent scope

    - localvars: VarMap
        Scope-local variable map

    - loc: Loc
        Start of scope location

    """

    def __init__(self, parent, loc):
        assert parent is None or isinstance(parent, Scope)
        assert isinstance(loc, Loc)
        self.parent = parent
        self.localvars = VarMap()
        self.loc = loc
        self.redefined = defaultdict(int)
        self.var_redefinitions = defaultdict(set)

    def define(self, name, loc):
        """
        Define a variable
        """
        v = Var(scope=self, name=name, loc=loc)
        self.localvars.define(v.name, v)
        return v

    def get(self, name):
        """
        Refer to a variable.  Returns the latest version.
        """
        if name in self.redefined:
            name = "%s.%d" % (name, self.redefined[name])
        return self.get_exact(name)

    def get_exact(self, name):
        """
        Refer to a variable.  The returned variable has the exact
        name (exact variable version).
        """
        try:
            return self.localvars.get(name)
        except NotDefinedError:
            if self.has_parent:
                return self.parent.get(name)
            else:
                raise

    def get_or_define(self, name, loc):
        if name in self.redefined:
            name = "%s.%d" % (name, self.redefined[name])

        if name not in self.localvars:
            return self.define(name, loc)
        else:
            return self.localvars.get(name)

    def redefine(self, name, loc, rename=True):
        """
        Redefine if the name is already defined
        """
        if name not in self.localvars:
            return self.define(name, loc)
        elif not rename:
            # Must use the same name if the variable is a cellvar, which
            # means it could be captured in a closure.
            return self.localvars.get(name)
        else:
            while True:
                ct = self.redefined[name]
                self.redefined[name] = ct + 1
                newname = "%s.%d" % (name, ct + 1)
                try:
                    res = self.define(newname, loc)
                except RedefinedError:
                    continue
                else:
                    self.var_redefinitions[name].add(newname)
                return res

    def get_versions_of(self, name):
        """
        Gets all known versions of a given name
        """
        vers = set()
        def walk(thename):
            redefs = self.var_redefinitions.get(thename, None)
            if redefs:
                for v in redefs:
                    vers.add(v)
                    walk(v)
        walk(name)
        return vers

    def make_temp(self, loc):
        n = len(self.localvars)
        v = Var(scope=self, name='$%d' % n, loc=loc)
        self.localvars.define(v.name, v)
        return v

    @property
    def has_parent(self):
        return self.parent is not None

    def __repr__(self):
        return "Scope(has_parent=%r, num_vars=%d, %s)" % (self.has_parent,
                                                          len(self.localvars),
                                                          self.loc)


class Block(EqualityCheckMixin):
    """A code block

    """

    def __init__(self, scope, loc):
        assert isinstance(scope, Scope)
        assert isinstance(loc, Loc)
        self.scope = scope
        self.body = []
        self.loc = loc

    def copy(self):
        block = Block(self.scope, self.loc)
        block.body = self.body[:]
        return block

    def find_exprs(self, op=None):
        """
        Iterate over exprs of the given *op* in this block.
        """
        for inst in self.body:
            if isinstance(inst, Assign):
                expr = inst.value
                if isinstance(expr, Expr):
                    if op is None or expr.op == op:
                        yield expr

    def find_insts(self, cls=None):
        """
        Iterate over insts of the given class in this block.
        """
        for inst in self.body:
            if isinstance(inst, cls):
                yield inst

    def find_variable_assignment(self, name):
        """
        Returns the assignment inst associated with variable "name", None if
        it cannot be found.
        """
        for x in self.find_insts(cls=Assign):
            if x.target.name == name:
                return x
        return None

    def prepend(self, inst):
        assert isinstance(inst, Stmt)
        self.body.insert(0, inst)

    def append(self, inst):
        assert isinstance(inst, Stmt)
        self.body.append(inst)

    def remove(self, inst):
        assert isinstance(inst, Stmt)
        del self.body[self.body.index(inst)]

    def clear(self):
        del self.body[:]

    def dump(self, file=None):
        # Avoid early bind of sys.stdout as default value
        file = file or sys.stdout
        for inst in self.body:
            if hasattr(inst, 'dump'):
                inst.dump(file)
            else:
                inst_vars = sorted(str(v) for v in inst.list_vars())
                print('    %-40s %s' % (inst, inst_vars), file=file)

    @property
    def terminator(self):
        return self.body[-1]

    @property
    def is_terminated(self):
        return self.body and self.body[-1].is_terminator

    def verify(self):
        if not self.is_terminated:
            raise VerificationError("Missing block terminator")
            # Only the last instruction can be a terminator
        for inst in self.body[:-1]:
            if inst.is_terminator:
                raise VerificationError("Terminator before the last "
                                        "instruction")

    def insert_after(self, stmt, other):
        """
        Insert *stmt* after *other*.
        """
        index = self.body.index(other)
        self.body.insert(index + 1, stmt)

    def insert_before_terminator(self, stmt):
        assert isinstance(stmt, Stmt)
        assert self.is_terminated
        self.body.insert(-1, stmt)

    def __repr__(self):
        return "<ir.Block at %s>" % (self.loc,)


class Loop(SlotEqualityCheckMixin):
    """Describes a loop-block
    """
    __slots__ = "entry", "exit"

    def __init__(self, entry, exit):
        self.entry = entry
        self.exit = exit

    def __repr__(self):
        args = self.entry, self.exit
        return "Loop(entry=%s, exit=%s)" % args


class With(SlotEqualityCheckMixin):
    """Describes a with-block
    """
    __slots__ = "entry", "exit"

    def __init__(self, entry, exit):
        self.entry = entry
        self.exit = exit

    def __repr__(self):
        args = self.entry, self.exit
        return "With(entry=%s, exit=%s)" % args


class FunctionIR(object):

    def __init__(self, blocks, is_generator, func_id, loc,
                 definitions, arg_count, arg_names):
        self.blocks = blocks
        self.is_generator = is_generator
        self.func_id = func_id
        self.loc = loc
        self.arg_count = arg_count
        self.arg_names = arg_names

        self._definitions = definitions

        self._reset_analysis_variables()

    def equal_ir(self, other):
        """ Checks that the IR contained within is equal to the IR in other.
        Equality is defined by being equal in fundamental structure (blocks,
        labels, IR node type and the order in which they are defined) and the
        IR nodes being equal. IR node equality essentially comes down to
        ensuring a node's `.__dict__` or `.__slots__` is equal, with the
        exception of ignoring 'loc' and 'scope' entries. The upshot is that the
        comparison is essentially location and scope invariant, but otherwise
        behaves as unsurprisingly as possible.
        """
        if type(self) is type(other):
            return self.blocks == other.blocks
        return False

    def diff_str(self, other):
        """
        Compute a human readable difference in the IR, returns a formatted
        string ready for printing.
        """
        msg = []
        for label, block in self.blocks.items():
            other_blk = other.blocks.get(label, None)
            if other_blk is not None:
                if block != other_blk:
                    msg.append(("Block %s differs" % label).center(80, '-'))
                    # see if the instructions are just a permutation
                    block_del = [x for x in block.body if isinstance(x, Del)]
                    oth_del = [x for x in other_blk.body if isinstance(x, Del)]
                    if block_del != oth_del:
                        # this is a common issue, dels are all present, but
                        # order shuffled.
                        if sorted(block_del) == sorted(oth_del):
                            msg.append(("Block %s contains the same dels but "
                                        "their order is different") % label)
                    if len(block.body) > len(other_blk.body):
                        msg.append("This block contains more statements")
                    elif len(block.body) < len(other_blk.body):
                        msg.append("Other block contains more statements")

                    # find the indexes where they don't match
                    tmp = []
                    for idx, stmts in enumerate(zip(block.body,
                                                    other_blk.body)):
                        b_s, o_s = stmts
                        if b_s != o_s:
                            tmp.append(idx)

                    def get_pad(ablock, l):
                        pointer = '-> '
                        sp = len(pointer) * ' '
                        pad = []
                        nstmt = len(ablock)
                        for i in range(nstmt):
                            if i in tmp:
                                item = pointer
                            elif i >= l:
                                item = pointer
                            else:
                                item = sp
                            pad.append(item)
                        return pad

                    min_stmt_len = min(len(block.body), len(other_blk.body))

                    with StringIO() as buf:
                        it = [("self", block), ("other", other_blk)]
                        for name, _block in it:
                            buf.truncate(0)
                            _block.dump(file=buf)
                            stmts = buf.getvalue().splitlines()
                            pad = get_pad(_block.body, min_stmt_len)
                            title = ("%s: block %s" % (name, label))
                            msg.append(title.center(80, '-'))
                            msg.extend(["{0}{1}".format(a, b) for a, b in
                                        zip(pad, stmts)])
        if msg == []:
            msg.append("IR is considered equivalent.")
        return '\n'.join(msg)

    def _reset_analysis_variables(self):

        self._consts = consts.ConstantInference(self)

        # Will be computed by PostProcessor
        self.generator_info = None
        self.variable_lifetime = None
        # { ir.Block: { variable names (potentially) alive at start of block } }
        self.block_entry_vars = {}

    def derive(self, blocks, arg_count=None, arg_names=None,
               force_non_generator=False, loc=None):
        """
        Derive a new function IR from this one, using the given blocks,
        and possibly modifying the argument count and generator flag.

        Post-processing will have to be run again on the new IR.
        """
        firstblock = blocks[min(blocks)]

        new_ir = copy.copy(self)
        new_ir.blocks = blocks
        new_ir.loc = firstblock.loc if loc is None else loc
        if force_non_generator:
            new_ir.is_generator = False
        if arg_count is not None:
            new_ir.arg_count = arg_count
        if arg_names is not None:
            new_ir.arg_names = arg_names
        new_ir._reset_analysis_variables()
        # Make fresh func_id
        new_ir.func_id = new_ir.func_id.derive()
        return new_ir

    def copy(self):
        new_ir = copy.copy(self)
        blocks = {}
        block_entry_vars = {}
        for label, block in self.blocks.items():
            new_block = block.copy()
            blocks[label] = new_block
            if block in self.block_entry_vars:
                block_entry_vars[new_block] = self.block_entry_vars[block]
        new_ir.blocks = blocks
        new_ir.block_entry_vars = block_entry_vars
        return new_ir

    def get_block_entry_vars(self, block):
        """
        Return a set of variable names possibly alive at the beginning of
        the block.
        """
        return self.block_entry_vars[block]

    def infer_constant(self, name):
        """
        Try to infer the constant value of a given variable.
        """
        if isinstance(name, Var):
            name = name.name
        return self._consts.infer_constant(name)

    def get_definition(self, value, lhs_only=False):
        """
        Get the definition site for the given variable name or instance.
        A Expr instance is returned by default, but if lhs_only is set
        to True, the left-hand-side variable is returned instead.
        """
        lhs = value
        while True:
            if isinstance(value, Var):
                lhs = value
                name = value.name
            elif isinstance(value, str):
                lhs = value
                name = value
            else:
                return lhs if lhs_only else value
            defs = self._definitions[name]
            if len(defs) == 0:
                raise KeyError("no definition for %r"
                               % (name,))
            if len(defs) > 1:
                raise KeyError("more than one definition for %r"
                               % (name,))
            value = defs[0]

    def get_assignee(self, rhs_value, in_blocks=None):
        """
        Finds the assignee for a given RHS value. If in_blocks is given the
        search will be limited to the specified blocks.
        """
        if in_blocks is None:
            blocks = self.blocks.values()
        elif isinstance(in_blocks, int):
            blocks = [self.blocks[in_blocks]]
        else:
            blocks = [self.blocks[blk] for blk in list(in_blocks)]

        assert isinstance(rhs_value, AbstractRHS)

        for blk in blocks:
            for assign in blk.find_insts(Assign):
                if assign.value == rhs_value:
                    return assign.target

        raise ValueError("Could not find an assignee for %s" % rhs_value)


    def dump(self, file=None):
        nofile = file is None
        # Avoid early bind of sys.stdout as default value
        file = file or StringIO()
        for offset, block in sorted(self.blocks.items()):
            print('label %s:' % (offset,), file=file)
            block.dump(file=file)
        if nofile:
            text = file.getvalue()
            if config.HIGHLIGHT_DUMPS:
                try:
                    import pygments
                except ImportError:
                    msg = "Please install pygments to see highlighted dumps"
                    raise ValueError(msg)
                else:
                    from pygments import highlight
                    from numba.misc.dump_style import NumbaIRLexer as lexer
                    from numba.misc.dump_style import by_colorscheme
                    from pygments.formatters import Terminal256Formatter
                    print(highlight(text, lexer(), Terminal256Formatter(
                        style=by_colorscheme())))
            else:
                print(text)


    def dump_to_string(self):
        with StringIO() as sb:
            self.dump(file=sb)
            return sb.getvalue()

    def dump_generator_info(self, file=None):
        file = file or sys.stdout
        gi = self.generator_info
        print("generator state variables:", sorted(gi.state_vars), file=file)
        for index, yp in sorted(gi.yield_points.items()):
            print("yield point #%d: live variables = %s, weak live variables = %s"
                  % (index, sorted(yp.live_vars), sorted(yp.weak_live_vars)),
                  file=file)

    def render_dot(self, filename_prefix="numba_ir", include_ir=True):
        """Render the CFG of the IR with GraphViz DOT via the
        ``graphviz`` python binding.

        Returns
        -------
        g : graphviz.Digraph
            Use `g.view()` to open the graph in the default PDF application.
        """

        try:
            import graphviz as gv
        except ImportError:
            raise ImportError(
                "The feature requires `graphviz` but it is not available. "
                "Please install with `pip install graphviz`"
            )
        g = gv.Digraph(
            filename="{}{}.dot".format(
                filename_prefix,
                self.func_id.unique_name,
            )
        )
        # Populate the nodes
        for k, blk in self.blocks.items():
            with StringIO() as sb:
                blk.dump(sb)
                label = sb.getvalue()
            if include_ir:
                label = ''.join(
                    [r'  {}\l'.format(x) for x in label.splitlines()],
                )
                label = r"block {}\l".format(k) + label
                g.node(str(k), label=label, shape='rect')
            else:
                label = r"{}\l".format(k)
                g.node(str(k), label=label, shape='circle')
        # Populate the edges
        for src, blk in self.blocks.items():
            for dst in blk.terminator.get_targets():
                g.edge(str(src), str(dst))
        return g


# A stub for undefined global reference
class UndefinedType(EqualityCheckMixin):

    _singleton = None

    def __new__(cls):
        obj = cls._singleton
        if obj is not None:
            return obj
        else:
            obj = object.__new__(cls)
            cls._singleton = obj
        return obj

    def __repr__(self):
        return "Undefined"


UNDEFINED = UndefinedType()
