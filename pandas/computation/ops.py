import operator as op

import numpy as np
from pandas.util.py3compat import PY3
import pandas.core.common as com
from pandas.core.base import StringMixin
from pandas.computation.common import flatten


_reductions = 'sum', 'prod'
_mathops = ('sin', 'cos', 'exp', 'log', 'expm1', 'log1p', 'pow', 'div', 'sqrt',
            'inv', 'sinh', 'cosh', 'tanh', 'arcsin', 'arccos', 'arctan',
            'arccosh', 'arcsinh', 'arctanh', 'arctan2', 'abs')


class OperatorError(Exception):
    pass


class UnaryOperatorError(OperatorError):
    pass


class BinaryOperatorError(OperatorError):
    pass


class Term(StringMixin):
    def __init__(self, name, env, side=None):
        self.name = name
        self.env = env
        self.side = side
        self.value = self._resolve_name()

        try:
            # ndframe potentially very slow for large, mixed dtype frames
            self.type = self.value.values.dtype
        except AttributeError:
            try:
                # ndarray
                self.type = self.value.dtype
            except AttributeError:
                # scalar
                self.type = type(self.value)

    def __unicode__(self):
        return com.pprint_thing(self.name)

    def _resolve_name(self):
        env = self.env
        key = self.name
        res = env.locals.get(key, env.globals.get(key))

        if res is None:
            if not isinstance(key, basestring):
                return key

            raise NameError('name {0!r} is not defined'.format(key))
        return res

    def update(self, value):
        env = self.env
        key = self.name
        if isinstance(key, basestring):
            try:
                del env.locals[key]
                env.locals[key] = value
            except KeyError:
                try:
                    del env.globals[key]
                    env.globals[key] = value
                except KeyError:
                    raise NameError('{0!r} is undefined'.format(key))

        self.value = value

    @property
    def isscalar(self):
        return np.isscalar(self.value)


class Constant(Term):
    def __init__(self, value, env):
        super(Constant, self).__init__(value, env)


def _print_operand(opr):
    return opr.name if is_term(opr) else unicode(opr)


class Op(StringMixin):
    """Hold an operator of unknown arity
    """
    def __init__(self, op, operands):
        self.op = op
        self.operands = operands

    def __iter__(self):
        return iter(self.operands)

    def __unicode__(self):
        """Print a generic n-ary operator and its operands using infix
        notation"""
        # recurse over the operands
        parened = ('({0})'.format(_print_operand(opr))
                   for opr in self.operands)
        return com.pprint_thing(' {0} '.format(self.op).join(parened))

    @property
    def return_type(self):
        # clobber types to bool if the op is a boolean operator
        if self.op in (_cmp_ops_syms + _bool_ops_syms):
            return np.bool_
        return np.result_type(*(term.type for term in flatten(self)))


_cmp_ops_syms = '>', '<', '>=', '<=', '==', '!=', '='
_cmp_ops_funcs = op.gt, op.lt, op.ge, op.le, op.eq, op.ne, op.eq
_cmp_ops_dict = dict(zip(_cmp_ops_syms, _cmp_ops_funcs))

_bool_ops_syms = '&', '|'
_bool_ops_funcs = op.and_, op.or_
_bool_ops_dict = dict(zip(_bool_ops_syms, _bool_ops_funcs))

_arith_ops_syms = '+', '-', '*', '/', '**', '//', '%'
_arith_ops_funcs = (op.add, op.sub, op.mul, op.truediv if PY3 else op.div,
                    op.pow, op.floordiv, op.mod)
_arith_ops_dict = dict(zip(_arith_ops_syms, _arith_ops_funcs))

_binary_ops_dict = {}

for d in (_cmp_ops_dict, _bool_ops_dict, _arith_ops_dict):
    _binary_ops_dict.update(d)


def _cast_inplace(terms, dtype):
    dt = np.dtype(dtype)
    for term in terms:
        # cast all the way down the tree since operands must be
        try:
            _cast_inplace(term.operands, dtype)
        except AttributeError:
            # we've bottomed out so actually do the cast
            try:
                new_value = term.value.astype(dt)
            except AttributeError:
                new_value = dt.type(term.value)
            term.update(new_value)


def is_term(obj):
    return isinstance(obj, Term)


def is_const(obj):
    return isinstance(obj, Constant)


class BinOp(Op):
    """Hold a binary operator and its operands

    Parameters
    ----------
    op : str or Op
    left : str or Op
    right : str or Op
    """
    def __init__(self, op, lhs, rhs):
        super(BinOp, self).__init__(op, (lhs, rhs))
        self.lhs = lhs
        self.rhs = rhs

        try:
            self.func = _binary_ops_dict[op]
        except KeyError:
            keys = _binary_ops_dict.keys()
            raise BinaryOperatorError('Invalid binary operator {0}, valid'
                                      ' operators are {1}'.format(op, keys))

    def __call__(self, env):
        # handle truediv
        if self.op == '/' and env.locals['truediv']:
            self.func = op.truediv

        # recurse over the left nodes
        try:
            left = self.lhs(env)
        except TypeError:
            left = self.lhs

        # recurse over the right nodes
        try:
            right = self.rhs(env)
        except TypeError:
            right = self.rhs

        # base cases
        if is_term(left) and is_term(right):
            res = self.func(left.value, right.value)
        elif not is_term(left) and is_term(right):
            res = self.func(left, right.value)
        elif is_term(left) and not is_term(right):
            res = self.func(left.value, right)
        elif not (is_term(left) or is_term(right)):
            res = self.func(left, right)

        return res


class Mod(BinOp):
    def __init__(self, lhs, rhs):
        super(Mod, self).__init__('%', lhs, rhs)
        _cast_inplace(self.operands, np.float_)


_unary_ops_syms = '+', '-', '~'
_unary_ops_funcs = op.pos, op.neg, op.invert
_unary_ops_dict = dict(zip(_unary_ops_syms, _unary_ops_funcs))


class UnaryOp(Op):
    """Hold a unary operator and its operands
    """
    def __init__(self, op, operand):
        super(UnaryOp, self).__init__(op, (operand,))
        self.operand = operand

        try:
            self.func = _unary_ops_dict[op]
        except KeyError:
            raise UnaryOperatorError('Invalid unary operator {0}, valid '
                                     'operators are '
                                     '{1}'.format(op, _unary_ops_syms))

    def __call__(self, env):
        operand = self.operand

        # recurse if operand is an Op
        try:
            operand = self.operand(env)
        except TypeError:
            operand = self.operand

        v = operand.value if is_term(operand) else operand

        try:
            res = self.func(v)
        except TypeError:
            res = self.func(v.values)

        return res

    def __unicode__(self):
        return com.pprint_thing('{0}({1})'.format(self.op, self.operand))

