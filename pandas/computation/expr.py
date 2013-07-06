import ast
import operator
import sys
import inspect
import itertools
import tokenize
from cStringIO import StringIO
from functools import partial

from pandas.core.base import StringMixin
from pandas.computation.ops import BinOp, UnaryOp, _reductions, _mathops
from pandas.computation.ops import _cmp_ops_syms, _bool_ops_syms
from pandas.computation.ops import _arith_ops_syms, _unary_ops_syms
from pandas.computation.ops import Term, Constant

import pandas.lib as lib
import datetime


class Scope(object):
    __slots__ = ('globals', 'locals', 'resolvers', '_global_resolvers',
                 'resolver_keys', '_resolver')

    def __init__(self, gbls=None, lcls=None, frame_level=1, resolvers=None):
        frame = sys._getframe(frame_level)

        try:
            self.globals = gbls or frame.f_globals.copy()
            self.locals = lcls or frame.f_locals.copy()
        finally:
            del frame

        # add some useful defaults
        self.globals['Timestamp'] = lib.Timestamp
        self.globals['datetime'] = datetime

        self.resolvers = resolvers or []
        self.resolver_keys = set(reduce(operator.add, (list(o.keys()) for o in
                                                       self.resolvers), []))
        self._global_resolvers = self.resolvers + [self.locals, self.globals]
        self._resolver = None

    @property
    def resolver(self):
        if self._resolver is None:
            def resolve_key(key):
                for resolver in self._global_resolvers:
                    try:
                        return resolver[key]
                    except KeyError:
                        pass
            self._resolver = resolve_key

        return self._resolver

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
            while sl >= 0:
                frame = frame.f_back
                sl -= 1
                frames.append(frame)
            for f in frames[::-1]:
                self.locals.update(f.f_locals)
        finally:
            del frame
            del frames


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


def _preparse(source):
    return _parenthesize_booleans(_rewrite_assign(source))



# partition all AST nodes
_all_nodes = frozenset(filter(lambda x: isinstance(x, type) and
                              issubclass(x, ast.AST),
                              (getattr(ast, node) for node in dir(ast))))


def _filter_nodes(superclass, all_nodes=_all_nodes):
    node_names = (node.__name__ for node in all_nodes
                  if issubclass(node, superclass))
    return frozenset(node_names)


_all_node_names = frozenset(map(lambda x: x.__name__, _all_nodes))
_mod_nodes = _filter_nodes(ast.mod)
_stmt_nodes = _filter_nodes(ast.stmt)
_expr_nodes = _filter_nodes(ast.expr)
_expr_context_nodes = _filter_nodes(ast.expr_context)
_slice_nodes = _filter_nodes(ast.slice)
_boolop_nodes = _filter_nodes(ast.boolop)
_operator_nodes = _filter_nodes(ast.operator)
_unary_op_nodes = _filter_nodes(ast.unaryop)
_cmp_op_nodes = _filter_nodes(ast.cmpop)
_comprehension_nodes = _filter_nodes(ast.comprehension)
_handler_nodes = _filter_nodes(ast.excepthandler)
_arguments_nodes = _filter_nodes(ast.arguments)
_keyword_nodes = _filter_nodes(ast.keyword)
_alias_nodes = _filter_nodes(ast.alias)


# nodes that we don't support directly but are needed for parsing
_hacked_nodes = frozenset(['Assign', 'Module', 'Expr'])


# these nodes are low priority or won't ever be supported (e.g., AST)
_unsupported_nodes = ((_stmt_nodes | _mod_nodes | _handler_nodes |
                       _arguments_nodes | _keyword_nodes | _alias_nodes |
                       _expr_context_nodes | frozenset(['Yield',
                                                        'GeneratorExp',
                                                        'IfExp', 'DictComp',
                                                        'SetComp', 'Repr',
                                                        'Lambda', 'Set', 'In',
                                                        'NotIn', 'AST',
                                                        'Is', 'IsNot'])) -
                      _hacked_nodes)

# we're adding a different assignment in some cases to be equality comparison
# and we don't want `stmt` and friends in their so get only the class whose
# names are capitalized
_base_supported_nodes = (_all_node_names - _unsupported_nodes) | _hacked_nodes
_msg = 'cannot both support and not support {0}'.format(_unsupported_nodes &
                                                        _base_supported_nodes)
assert not _unsupported_nodes & _base_supported_nodes, _msg


def _node_not_implemented(node_name, cls):
    def f(self, *args, **kwargs):
        raise NotImplementedError("{0!r} nodes are not "
                                  "implemented".format(node_name))
    return f


def disallow(nodes):
    def disallowed(cls):
        cls.unsupported_nodes = ()
        for node in nodes:
            new_method =  _node_not_implemented(node, cls)
            name = 'visit_{0}'.format(node)
            cls.unsupported_nodes += (name,)
            setattr(cls, name, new_method)
        return cls
    return disallowed


def _op_maker(op_class, op_symbol):
    def f(self, node, *args, **kwargs):
        return partial(op_class, op_symbol, *args, **kwargs)
    return f


_op_classes = {'binary': BinOp, 'unary': UnaryOp}

def add_ops(op_classes):
    def f(cls):
        for op_attr_name, op_class in op_classes.iteritems():
            ops = getattr(cls, '{0}_ops'.format(op_attr_name))
            ops_map = getattr(cls, '{0}_op_nodes_map'.format(op_attr_name))
            for op in ops:
                setattr(cls, 'visit_{0}'.format(ops_map[op]),
                        _op_maker(op_class, op))
        return cls
    return f


@disallow(_unsupported_nodes)
@add_ops(_op_classes)
class BaseExprVisitor(ast.NodeVisitor):

    """Custom ast walker
    """
    binary_ops = _cmp_ops_syms + _bool_ops_syms + _arith_ops_syms
    binary_op_nodes = ('Gt', 'Lt', 'GtE', 'LtE', 'Eq', 'NotEq', 'BitAnd',
                       'BitOr', 'Add', 'Sub', 'Mult', 'Div', 'Pow', 'FloorDiv',
                       'Mod')
    binary_op_nodes_map = dict(itertools.izip(binary_ops, binary_op_nodes))

    unary_ops = _unary_ops_syms
    unary_op_nodes = 'UAdd', 'USub', 'Invert'
    unary_op_nodes_map = dict(itertools.izip(unary_ops, unary_op_nodes))

    def __init__(self, env, preparser=_preparse):
        self.env = env
        self.preparser = preparser

    def visit(self, node, **kwargs):
        if isinstance(node, basestring):
            node = ast.fix_missing_locations(ast.parse(self.preparser(node)))

        method = 'visit_' + node.__class__.__name__
        visitor = getattr(self, method, None)
        return visitor(node, **kwargs)

    def visit_Module(self, node, **kwargs):
        if len(node.body) != 1:
            raise SyntaxError('only a single expression is allowed')
        expr = node.body[0]
        return self.visit(expr, **kwargs)

    def visit_Expr(self, node, **kwargs):
        return self.visit(node.value, **kwargs)

    def visit_BinOp(self, node, **kwargs):
        op = self.visit(node.op)
        left = self.visit(node.left, side='left')
        right = self.visit(node.right, side='right')
        return op(left, right)

    def visit_UnaryOp(self, node, **kwargs):
        op = self.visit(node.op)
        return op(self.visit(node.operand))

    def visit_Name(self, node, **kwargs):
        return Term(node.id, self.env)

    def visit_Num(self, node, **kwargs):
        return Constant(node.n, self.env)

    def visit_Str(self, node, **kwargs):
        return Constant(node.s, self.env)

    def visit_List(self, node, **kwargs):
        return Constant([self.visit(e).value for e in node.elts], self.env)

    visit_Tuple = visit_List

    def visit_Index(self, node, **kwargs):
        """ df.index[4] """
        return self.visit(node.value)


    def visit_Subscript(self, node, **kwargs):
        """ df.index[4:6] """
        value = self.visit(node.value)
        slobj = self.visit(node.slice)

        try:
            return Constant(value[slobj], self.env)
        except TypeError:
            raise ValueError("cannot subscript [{0}] with "
                             "[{1}]".format(value, slobj))

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

        return slice(lower, upper, step)

    def visit_Assign(self, node, **kwargs):
        cmpr = ast.Compare(ops=[ast.Eq()], left=node.targets[0],
                           comparators=[node.value])
        return self.visit(cmpr)

    def visit_Attribute(self, node, **kwargs):
        attr = node.attr
        value = node.value

        ctx = node.ctx.__class__
        if ctx == ast.Load:
            # resolve the value
            return getattr(self.visit(value).value, attr)
        raise ValueError("Invalid Attribute context {0}".format(ctx.__name__))

    def visit_Call(self, node, **kwargs):

        # this can happen with: datetime.datetime
        if isinstance(node.func, ast.Attribute):
            res = self.visit_Attribute(node.func)
        elif not isinstance(node.func, ast.Name):
            raise TypeError("Only named functions are supported")
        else:
            res = self.visit(node.func)

        if res is None:
            raise ValueError("Invalid function call {0}".format(node.func.id))
        if hasattr(res, 'value'):
            res = res.value

        args = [self.visit(targ).value for targ in node.args]
        if node.starargs is not None:
            args = args + self.visit(node.starargs).value

        keywords = {}
        for key in node.keywords:
            if not isinstance(key, ast.keyword):
                raise ValueError(
                    "keyword error in function call '{0}'".format(node.func.id))
            keywords[key.arg] = self.visit(key.value).value
        if node.kwargs is not None:
            keywords.update(self.visit(node.kwargs).value)

        return Constant(res(*args, **keywords), self.env)

    def visit_Compare(self, node, **kwargs):
        ops = node.ops
        comps = node.comparators
        for op, comp in itertools.izip(ops, comps):
            vop = self.visit(op)
            node = vop(self.visit(node.left, side='left'),
                       self.visit(comp, side='right'))
        return node


_numexpr_not_supported = frozenset(['Assign', 'BoolOp', 'Not', 'Str', 'Slice',
                                    'Index', 'Subscript', 'Tuple', 'List',
                                    'Dict', 'Call'])
_numexpr_supported_calls = frozenset(_reductions + _mathops)

@disallow(_unsupported_nodes | _numexpr_not_supported)
class NumExprVisitor(BaseExprVisitor):
    def __init__(self, env, preparser=None):
        if preparser is not None:
            raise ValueError("only strict numexpr syntax is supported")
        preparser = lambda x: x
        super(NumExprVisitor, self).__init__(env, preparser)


_python_not_supported = _numexpr_not_supported

@disallow(_unsupported_nodes | _python_not_supported)
class PythonExprVisitor(BaseExprVisitor):
    pass


class Expr(StringMixin):

    """Expr object"""

    def __init__(self, expr, engine='numexpr', env=None, truediv=True):
        self.expr = expr
        self.env = env or Scope(frame_level=2)
        self._visitor = _visitors[engine](self.env)
        self.terms = self.parse()
        self.engine = engine
        self.truediv = truediv

    def __call__(self, env):
        env.locals['truediv'] = self.truediv
        return self.terms(env)

    def __unicode__(self):
        return unicode(self.terms)

    def __len__(self):
        return len(self.expr)

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


_visitors = {'python': PythonExprVisitor, 'numexpr': NumExprVisitor}
