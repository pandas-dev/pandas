import operator as op
from functools import partial

from pandas.util.py3compat import PY3


_reductions = 'sum', 'prod'
_mathops = 'sin', 'cos', 'tan'


class OperatorError(Exception):
    pass


class UnaryOperatorError(OperatorError):
    pass


class BinaryOperatorError(OperatorError):
    pass


def _resolve_name(env, key):
    res = env.locals.get(key, env.globals.get(key))

    if res is None:
        if not isinstance(key, basestring):
            return key

        raise NameError('{0!r} is undefined'.format(key))

    return res


def _update_name(env, key, value):
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


def _update_names(env, mapping):
    updater = partial(_update_name, env)
    for key, value in mapping.iteritems():
        updater(key, value)


class Op(object):
    """Hold an operator of unknown arity
    """
    def __init__(self, op, operands):
        self.op = op
        self.operands = operands

    def __iter__(self):
        return iter(self.operands)

    @property
    def name(self):
        return self.__class__.__name__


_cmp_ops_syms = '>', '<', '>=', '<=', '==', '!='
_cmp_ops_funcs = op.gt, op.lt, op.ge, op.le, op.eq, op.ne
_cmp_ops_dict = dict(zip(_cmp_ops_syms, _cmp_ops_funcs))

_bool_ops_syms = '&', '|'
_bool_ops_funcs = op.and_, op.or_
_bool_ops_dict = dict(zip(_bool_ops_syms, _bool_ops_funcs))

_arith_ops_syms = '+', '-', '*', '/', '**', '//'
_arith_ops_funcs = (op.add, op.sub, op.mul, op.truediv if PY3 else op.div,
                    op.pow, op.floordiv)
_arith_ops_dict = dict(zip(_arith_ops_syms, _arith_ops_funcs))

_binary_ops_dict = {}

for d in (_cmp_ops_dict, _bool_ops_dict, _arith_ops_dict):
    _binary_ops_dict.update(d)


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

    def __repr__(self):
        return '{0}(op={1!r}, lhs={2!r}, rhs={3!r})'.format(self.name, self.op,
                                                            self.lhs, self.rhs)

    __str__ = __repr__

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
        if not (isinstance(left, basestring) or isinstance(right, basestring)):
            res = self.func(left, right)
        elif isinstance(left, basestring) and not isinstance(right,
                                                             basestring):
            res = self.func(_resolve_name(env, left), right)
        elif not isinstance(left, basestring) and isinstance(right,
                                                             basestring):
            res = self.func(left, _resolve_name(env, right))
        elif isinstance(left, basestring) and isinstance(right, basestring):
            res = self.func(_resolve_name(env, left), _resolve_name(env,
                                                                    right))

        return res


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

        if isinstance(operand, basestring):
            v = _resolve_name(env, operand)
        else:
            v = operand

        try:
            res = self.func(v)
        except TypeError:
            res = self.func(v.values)

        return res

    def __repr__(self):
        return '{0}(op={1!r}, operand={2!r})'.format(self.name, self.op,
                                                     self.operand)
