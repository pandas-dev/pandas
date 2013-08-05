import re
import operator as op
from functools import partial
from itertools import product, islice, chain

import numpy as np

import pandas as pd
from pandas.compat import PY3, string_types, text_type
import pandas.core.common as com
from pandas.core.base import StringMixin
from pandas.computation.common import _ensure_decoded


_reductions = 'sum', 'prod'
_mathops = ('sin', 'cos', 'exp', 'log', 'expm1', 'log1p', 'pow', 'div', 'sqrt',
            'inv', 'sinh', 'cosh', 'tanh', 'arcsin', 'arccos', 'arctan',
            'arccosh', 'arcsinh', 'arctanh', 'arctan2', 'abs')


_LOCAL_TAG = '__pd_eval_local_'
_TAG_RE = re.compile('^{0}'.format(_LOCAL_TAG))


class UndefinedVariableError(NameError):
    def __init__(self, *args):
        msg = 'name {0!r} is not defined'
        subbed = _TAG_RE.sub('', args[0])
        if subbed != args[0]:
            subbed = '@' + subbed
            msg = 'local variable {0!r} is not defined'
        super(UndefinedVariableError, self).__init__(msg.format(subbed))


class OperatorError(Exception):
    pass


class UnaryOperatorError(OperatorError):
    pass


class BinaryOperatorError(OperatorError):
    pass


def _possibly_update_key(d, value, old_key, new_key=None):
    if new_key is None:
        new_key = old_key

    try:
        del d[old_key]
    except KeyError:
        return False
    else:
        d[new_key] = value
        return True


class Term(StringMixin):
    def __init__(self, name, env, side=None, encoding=None):
        self._name = name
        self.env = env
        self.side = side
        self.local = _TAG_RE.search(text_type(name)) is not None
        self._value = self._resolve_name()
        self.encoding = encoding

    @property
    def local_name(self):
        return _TAG_RE.sub('', self.name)

    def __unicode__(self):
        return com.pprint_thing(self.name)

    def __call__(self, *args, **kwargs):
        return self.value

    def _resolve_name(self):
        #import ipdb; ipdb.set_trace()
        env = self.env
        key = self.name
        res = env.resolve(self.local_name, globally=not self.local)
        self.update(res)

        if res is None:
            if not isinstance(key, string_types):
                return key
            raise UndefinedVariableError(key)

        if hasattr(res, 'ndim') and res.ndim > 2:
            raise NotImplementedError("N-dimensional objects, where N > 2, are"
                                      " not supported with eval")
        return res

    def update(self, value):
        """
        search order for local (i.e., @variable) variables:

        scope, key_variable
        [('locals', 'local_name'),
         ('globals', 'local_name'),
         ('locals', 'key'),
         ('globals', 'key')]
        """
        env = self.env
        key = self.name

        # if it's a variable name (otherwise a constant)
        if isinstance(key, string_types):
            if self.local:
                # get it's name WITHOUT the local tag (defined above)
                local_name = self.local_name

                # search for the local in the above specified order
                scope_pairs = product([env.locals, env.globals],
                                      [local_name, key])

                # a[::2] + a[1::2] but iterators
                scope_iter = chain(islice(scope_pairs, None, None, 2),
                                   islice(scope_pairs, 1, None, 2))
                for d, k in scope_iter:
                    if _possibly_update_key(d, value, k, key):
                        break
                else:
                    raise UndefinedVariableError(key)
            else:
                # otherwise we look in resolvers -> locals -> globals
                for r in (env.resolver_dict, env.locals, env.globals):
                    if _possibly_update_key(r, value, key):
                        break
                else:
                    raise UndefinedVariableError(key)

        self.value = value

    @property
    def isscalar(self):
        return np.isscalar(self._value)

    @property
    def type(self):
        try:
            # potentially very slow for large, mixed dtype frames
            return self._value.values.dtype
        except AttributeError:
            try:
                # ndarray
                return self._value.dtype
            except AttributeError:
                # scalar
                return type(self._value)

    return_type = type

    @property
    def raw(self):
        return com.pprint_thing('{0}(name={1!r}, type={2})'
                                ''.format(self.__class__.__name__, self.name,
                                          self.type))

    @property
    def kind(self):
        try:
            return self.type.__name__
        except AttributeError:
            return self.type.type.__name__

    @property
    def value(self):
        kind = self.kind.lower()
        if kind == 'datetime64':
            try:
                return self._value.asi8
            except AttributeError:
                return self._value.view('i8')
        elif kind == 'datetime':
            return pd.Timestamp(self._value)
        elif kind == 'timestamp':
            return self._value.asm8.view('i8')
        return self._value

    @value.setter
    def value(self, new_value):
        self._value = new_value

    @property
    def name(self):
        return self._name

    @name.setter
    def name(self, new_name):
        self._name = new_name


class Constant(Term):
    def __init__(self, value, env):
        super(Constant, self).__init__(value, env)

    def _resolve_name(self):
        return self._name

    @property
    def name(self):
        return self.value


def _print_operand(opr):
    return opr.name if is_term(opr) else com.pprint_thing(opr)


def _get_op(op):
    return {'not': '~', 'and': '&', 'or': '|'}.get(op, op)


class Op(StringMixin):
    """Hold an operator of unknown arity
    """
    def __init__(self, op, operands, *args, **kwargs):
        self.op = _get_op(op)
        self.operands = operands
        self.encoding = kwargs.get('encoding', None)

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
        return np.result_type(*(term.type for term in com.flatten(self)))

    @property
    def raw(self):
        parened = ('{0}({1!r}, {2})'.format(self.__class__.__name__, self.op,
                                            ', '.join('{0}'.format(opr.raw) for
                                                      opr in self.operands)))
        return parened


_cmp_ops_syms = '>', '<', '>=', '<=', '==', '!='
_cmp_ops_funcs = op.gt, op.lt, op.ge, op.le, op.eq, op.ne
_cmp_ops_dict = dict(zip(_cmp_ops_syms, _cmp_ops_funcs))

_bool_ops_syms = '&', '|', 'and', 'or'
_bool_ops_funcs = op.and_, op.or_, op.and_, op.or_
_bool_ops_dict = dict(zip(_bool_ops_syms, _bool_ops_funcs))

_arith_ops_syms = '+', '-', '*', '/', '**', '//', '%'
_arith_ops_funcs = (op.add, op.sub, op.mul, op.truediv if PY3 else op.div,
                    op.pow, op.floordiv, op.mod)
_arith_ops_dict = dict(zip(_arith_ops_syms, _arith_ops_funcs))

_special_case_arith_ops_syms = '**', '//', '%'
_special_case_arith_ops_funcs = op.pow, op.floordiv, op.mod
_special_case_arith_ops_dict = dict(zip(_special_case_arith_ops_syms,
                                        _special_case_arith_ops_funcs))

_binary_ops_dict = {}

for d in (_cmp_ops_dict, _bool_ops_dict, _arith_ops_dict):
    _binary_ops_dict.update(d)


def _cast_inplace(terms, dtype):
    dt = np.dtype(dtype)
    for term in terms:
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
    def __init__(self, op, lhs, rhs, **kwargs):
        super(BinOp, self).__init__(op, (lhs, rhs))
        self.lhs = lhs
        self.rhs = rhs

        self.convert_values()

        try:
            self.func = _binary_ops_dict[op]
        except KeyError:
            keys = _binary_ops_dict.keys()
            raise BinaryOperatorError('Invalid binary operator {0!r}, valid'
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

    def convert_values(self):
        def stringify(value):
            if self.encoding is not None:
                encoder = partial(com.pprint_thing_encoded,
                                  encoding=self.encoding)
            else:
                encoder = com.pprint_thing
            return encoder(value)

        lhs, rhs = self.lhs, self.rhs

        if (is_term(lhs) and lhs.kind.startswith('datetime') and is_term(rhs)
                and rhs.isscalar):
            v = rhs.value
            if isinstance(v, (int, float)):
                v = stringify(v)
            v = _ensure_decoded(v)
            v = pd.Timestamp(v)
            if v.tz is not None:
                v = v.tz_convert('UTC')
            self.rhs.update(v)

        if (is_term(rhs) and rhs.kind.startswith('datetime') and
                is_term(lhs) and lhs.isscalar):
            v = lhs.value
            if isinstance(v, (int, float)):
                v = stringify(v)
            v = _ensure_decoded(v)
            v = pd.Timestamp(v)
            if v.tz is not None:
                v = v.tz_convert('UTC')
            self.lhs.update(v)


class Div(BinOp):
    def __init__(self, lhs, rhs, truediv=True, *args, **kwargs):
        super(Div, self).__init__('/', lhs, rhs, *args, **kwargs)
        if truediv or PY3:
            _cast_inplace(com.flatten(self), np.float_)


_unary_ops_syms = '+', '-', '~', 'not'
_unary_ops_funcs = op.pos, op.neg, op.invert, op.invert
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
