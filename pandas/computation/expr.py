import ast
import sys, inspect
import itertools
import tokenize
import re
from cStringIO import StringIO
from functools import partial

from pandas.core.base import StringMixin
from pandas.computation.ops import BinOp, UnaryOp, _reductions, _mathops
from pandas.computation.ops import _cmp_ops_syms, _bool_ops_syms
from pandas.computation.ops import _arith_ops_syms, _unary_ops_syms
from pandas.computation.ops import Term, Constant, Value

from pandas import Timestamp
import datetime

class Scope(object):
    __slots__ = 'globals', 'locals'

    def __init__(self, gbls=None, lcls=None, frame_level=1):
        frame = sys._getframe(frame_level)

        try:
            self.globals = gbls or frame.f_globals.copy()
            self.locals = lcls or frame.f_locals.copy()
        finally:
            del frame

        # add some usefule defaults
        self.globals['Timestamp'] = Timestamp
        self.globals['datetime'] = datetime

    def update(self, scope_level=None):

        # we are always 2 levels below the caller
        # plus the caller maybe below the env level
        # in which case we need addtl levels
        sl = 2
        if scope_level is not None:
            sl += scope_level

        # add sl frames to the scope starting with the
        # most distant and overwritting with more current
        # makes sure that we can capture variable scope
        frame = inspect.currentframe()
        try:
            frames = []
            while(sl>=0):
                frame = frame.f_back
                sl -= 1
                frames.append(frame)
            for f in frames[::-1]:
                self.locals.update(f.f_locals)
        finally:
            del frame
            del frames

class ExprParserError(Exception):
    pass


def _rewrite_assign(source):
    res = []
    g = tokenize.generate_tokens(StringIO(source).readline)
    for toknum, tokval, _, _, _ in g:
        res.append((toknum, '==' if tokval == '=' else tokval))
    return tokenize.untokenize(res)


def _parenthesize_booleans(source, ops='|&'):
    res = source
    for op in ops:
        terms = res.split(op)

        t = []
        for term in terms:
            t.append('({0})'.format(term))

        res = op.join(t)
    return res


def preparse(source):
    return _parenthesize_booleans(_rewrite_assign(source))


class ExprVisitor(ast.NodeVisitor):
    """Custom ast walker
    """
    bin_ops = _cmp_ops_syms + _bool_ops_syms + _arith_ops_syms
    bin_op_nodes = ('Gt', 'Lt', 'GtE', 'LtE', 'Eq', 'NotEq', None,
                    'BitAnd', 'BitOr', 'Add', 'Sub', 'Mult', 'Div', 'Pow',
                    'FloorDiv', 'Mod')
    bin_op_nodes_map = dict(zip(bin_ops, bin_op_nodes))

    unary_ops = _unary_ops_syms
    unary_op_nodes = 'UAdd', 'USub', 'Invert'
    unary_op_nodes_map = dict(zip(unary_ops, unary_op_nodes))

    def __init__(self, env):
        for bin_op in itertools.ifilter(lambda x: x is not None, self.bin_ops):
            setattr(self, 'visit_{0}'.format(self.bin_op_nodes_map[bin_op]),
                    lambda node, bin_op=bin_op: partial(BinOp, bin_op))

        for unary_op in self.unary_ops:
            setattr(self,
                    'visit_{0}'.format(self.unary_op_nodes_map[unary_op]),
                    lambda node, unary_op=unary_op: partial(UnaryOp, unary_op))
        self.env = env

    def not_implemented(self, s):
        raise NotImplementedError("{0} not yet supported".format(s))

    def visit(self, node, **kwargs):
        if not (isinstance(node, ast.AST) or isinstance(node, basestring)):
            raise TypeError('"node" must be an AST node or a string, you'
                            ' passed a(n) {0}'.format(node.__class__))
        if isinstance(node, basestring):
            node = ast.fix_missing_locations(ast.parse(preparse(node)))

        method = 'visit_' + node.__class__.__name__
        visitor = getattr(self, method, None)
        if visitor is None:
            self.not_implemented("ast visitor [{0}]".format(method))
        return visitor(node, **kwargs)

    def visit_Module(self, node, **kwargs):
        if len(node.body) != 1:
            raise ExprParserError('only a single expression is allowed')

        expr = node.body[0]
        if not isinstance(expr, (ast.Expr, ast.Assign)):
            raise SyntaxError('only expressions are allowed')

        return self.visit(expr, **kwargs)

    def visit_Expr(self, node, **kwargs):
        return self.visit(node.value, **kwargs)

    def visit_BinOp(self, node, **kwargs):
        op = self.visit(node.op)
        left = self.visit(node.left,side='left')
        right = self.visit(node.right,side='right')
        return op(left, right)

    def visit_UnaryOp(self, node, **kwargs):
        if isinstance(node.op, ast.Not):
            self.not_implemented('not operator')
        op = self.visit(node.op)
        return op(self.visit(node.operand))

    def visit_List(self, node, **kwargs):
        return Value([ self.visit(e).value for e in node.elts ], self.env)

    def visit_Name(self, node, **kwargs):
        return Term(node.id, self.env)

    def visit_Num(self, node, **kwargs):
        return Constant(node.n, self.env)

    def visit_Str(self, node, **kwargs):
        return Value(node.s, self.env)

    def visit_Index(self, node, **kwargs):
        """ df.index[4] """
        return self.visit(node.value).value

    def visit_Subscript(self, node, **kwargs):
        """ df.index[4:6] """
        value = self.visit(node.value)
        slobj = self.visit(node.slice)

        return Value(value[slobj],self.env)

    def visit_Slice(self, node, **kwargs):
        """ df.index[slice(4,6)] """
        lower = node.lower
        if lower is not None:
            lower = self.visit(lower).value
        upper = node.upper
        if upper is not None:
            upper = self.visit(upper).value
        step = node.step
        if step is not None:
            step = self.visit(step).value

        return slice(lower,upper,step)

    def visit_Compare(self, node, **kwargs):
        ops = node.ops
        comps = node.comparators
        if len(ops) != 1:
            raise ExprParserError('chained comparisons not supported')
        return self.visit(ops[0])(self.visit(node.left,side='left'), self.visit(comps[0],side='right'))

    def visit_Assign(self, node, **kwargs):
        cmpr = ast.copy_location(ast.Compare(ops=[ast.Eq()],
                                             left=node.targets[0],
                                             comparators=[node.value]), node)
        return self.visit(cmpr)

    def visit_Call(self, node, **kwargs):
        if not isinstance(node.func, ast.Name):
            raise TypeError("Only named functions are supported")

        valid_ops = _reductions + _mathops

        if node.func.id not in valid_ops:
            raise ValueError("Only {0} are supported".format(valid_ops))

        self.not_implemented('function calls')

    def visit_Attribute(self, node, **kwargs):
        self.not_implemented('attribute access')

    def visit_BoolOp(self, node, **kwargs):
        self.not_implemented('boolean operators')


class Expr(StringMixin):
    """Expr object"""
    def __init__(self, expr, engine='numexpr', env=None, truediv=True):
        self.expr = expr
        self.env = env or Scope(frame_level=2)
        self._visitor = ExprVisitor(self.env)
        self.terms = self.parse()
        self.engine = engine
        self.truediv = truediv

    def __call__(self, env):
        env.locals['truediv'] = self.truediv
        return self.terms(env)

    def __unicode__(self):
        return unicode(self.terms)

    def parse(self):
        """return a Termset"""
        return self._visitor.visit(self.expr)

    def align(self):
        """align a set of Terms"""
        return self.terms.align(self.env)


def isexpr(s, check_names=True):
    try:
        Expr(s)
    except SyntaxError:
        return False
    except NameError:
        return not check_names
    else:
        return True
