""":func:`~pandas.eval` parsers
"""

import ast
import operator
import sys
import inspect
import tokenize
import datetime
import struct

from functools import partial

import pandas as pd
from pandas import compat
from pandas.compat import StringIO, zip, reduce, string_types
from pandas.core.base import StringMixin
from pandas.core import common as com
from pandas.computation.common import NameResolutionError
from pandas.computation.ops import (_cmp_ops_syms, _bool_ops_syms,
                                    _arith_ops_syms, _unary_ops_syms, is_term)
from pandas.computation.ops import _reductions, _mathops, _LOCAL_TAG
from pandas.computation.ops import Op, BinOp, UnaryOp, Term, Constant, Div
from pandas.computation.ops import UndefinedVariableError


def _ensure_scope(level=2, global_dict=None, local_dict=None, resolvers=None,
                  target=None, **kwargs):
    """Ensure that we are grabbing the correct scope."""
    return Scope(gbls=global_dict, lcls=local_dict, level=level,
                 resolvers=resolvers, target=target)


def _check_disjoint_resolver_names(resolver_keys, local_keys, global_keys):
    """Make sure that variables in resolvers don't overlap with locals or
    globals.
    """
    res_locals = list(com.intersection(resolver_keys, local_keys))
    if res_locals:
        msg = "resolvers and locals overlap on names {0}".format(res_locals)
        raise NameResolutionError(msg)

    res_globals = list(com.intersection(resolver_keys, global_keys))
    if res_globals:
        msg = "resolvers and globals overlap on names {0}".format(res_globals)
        raise NameResolutionError(msg)


def _replacer(x, pad_size):
    """Replace a number with its padded hexadecimal representation. Used to tag
    temporary variables with their calling scope's id.
    """
    # get the hex repr of the binary char and remove 0x and pad by pad_size
    # zeros
    try:
        hexin = ord(x)
    except TypeError:
        # bytes literals masquerade as ints when iterating in py3
        hexin = x

    return hex(hexin).replace('0x', '').rjust(pad_size, '0')


def _raw_hex_id(obj, pad_size=2):
    """Return the padded hexadecimal id of ``obj``."""
    # interpret as a pointer since that's what really what id returns
    packed = struct.pack('@P', id(obj))

    return ''.join(_replacer(x, pad_size) for x in packed)


class Scope(StringMixin):

    """Object to hold scope, with a few bells to deal with some custom syntax
    added by pandas.

    Parameters
    ----------
    gbls : dict or None, optional, default None
    lcls : dict or Scope or None, optional, default None
    level : int, optional, default 1
    resolvers : list-like or None, optional, default None

    Attributes
    ----------
    globals : dict
    locals : dict
    level : int
    resolvers : tuple
    resolver_keys : frozenset
    """
    __slots__ = ('globals', 'locals', 'resolvers', '_global_resolvers',
                 'resolver_keys', '_resolver', 'level', 'ntemps', 'target')

    def __init__(self, gbls=None, lcls=None, level=1, resolvers=None,
                 target=None):
        self.level = level
        self.resolvers = tuple(resolvers or [])
        self.globals = dict()
        self.locals = dict()
        self.target = target
        self.ntemps = 1  # number of temporary variables in this scope

        if isinstance(lcls, Scope):
            ld, lcls = lcls, dict()
            self.locals.update(ld.locals.copy())
            self.globals.update(ld.globals.copy())
            self.resolvers += ld.resolvers
            if ld.target is not None:
                self.target = ld.target
            self.update(ld.level)

        frame = sys._getframe(level)
        try:
            self.globals.update(gbls or frame.f_globals)
            self.locals.update(lcls or frame.f_locals)
        finally:
            del frame

        # add some useful defaults
        self.globals['Timestamp'] = pd.lib.Timestamp
        self.globals['datetime'] = datetime

        # SUCH a hack
        self.globals['True'] = True
        self.globals['False'] = False

        # function defs
        self.globals['list'] = list
        self.globals['tuple'] = tuple

        res_keys = (list(o.keys()) for o in self.resolvers)
        self.resolver_keys = frozenset(reduce(operator.add, res_keys, []))
        self._global_resolvers = self.resolvers + (self.locals, self.globals)
        self._resolver = None

        self.resolver_dict = {}
        for o in self.resolvers:
            self.resolver_dict.update(dict(o))

    def __unicode__(self):
        return com.pprint_thing(
            'locals: {0}\nglobals: {0}\nresolvers: '
            '{0}\ntarget: {0}'.format(list(self.locals.keys()),
                                      list(self.globals.keys()),
                                      list(self.resolver_keys),
                                      self.target))

    def __getitem__(self, key):
        return self.resolve(key, globally=False)

    def resolve(self, key, globally=False):
        resolvers = self.locals, self.globals
        if globally:
            resolvers = self._global_resolvers

        for resolver in resolvers:
            try:
                return resolver[key]
            except KeyError:
                pass

    def update(self, level=None):
        """Update the current scope by going back `level` levels.

        Parameters
        ----------
        level : int or None, optional, default None
        """
        # we are always 2 levels below the caller
        # plus the caller may be below the env level
        # in which case we need addtl levels
        sl = 2
        if level is not None:
            sl += level

        # add sl frames to the scope starting with the
        # most distant and overwritting with more current
        # makes sure that we can capture variable scope
        frame = inspect.currentframe()
        try:
            frames = []
            while sl >= 0:
                frame = frame.f_back
                sl -= 1
                if frame is None:
                    break
                frames.append(frame)
            for f in frames[::-1]:
                self.locals.update(f.f_locals)
                self.globals.update(f.f_globals)
        finally:
            del frame, frames

    def add_tmp(self, value, where='locals'):
        """Add a temporary variable to the scope.

        Parameters
        ----------
        value : object
            An arbitrary object to be assigned to a temporary variable.
        where : basestring, optional, default 'locals', {'locals', 'globals'}
            What scope to add the value to.

        Returns
        -------
        name : basestring
            The name of the temporary variable created.
        """
        d = getattr(self, where, None)

        if d is None:
            raise AttributeError("Cannot add value to non-existent scope "
                                 "{0!r}".format(where))
        if not isinstance(d, dict):
            raise TypeError("Cannot add value to object of type {0!r}, "
                            "scope must be a dictionary"
                            "".format(type(d).__name__))
        name = 'tmp_var_{0}_{1}_{2}'.format(type(value).__name__, self.ntemps,
                                            _raw_hex_id(self))
        d[name] = value

        # only increment if the variable gets put in the scope
        self.ntemps += 1
        return name

    def remove_tmp(self, name, where='locals'):
        d = getattr(self, where, None)
        if d is None:
            raise AttributeError("Cannot remove value from non-existent scope "
                                 "{0!r}".format(where))
        if not isinstance(d, dict):
            raise TypeError("Cannot remove value from object of type {0!r}, "
                            "scope must be a dictionary"
                            "".format(type(d).__name__))
        del d[name]
        self.ntemps -= 1


def _rewrite_assign(source):
    """Rewrite the assignment operator for PyTables expression that want to use
    ``=`` as a substitute for ``==``.
    """
    res = []
    g = tokenize.generate_tokens(StringIO(source).readline)
    for toknum, tokval, _, _, _ in g:
        res.append((toknum, '==' if tokval == '=' else tokval))
    return tokenize.untokenize(res)


def _replace_booleans(source):
    """Replace ``&`` with ``and`` and ``|`` with ``or`` so that bitwise
    precedence is changed to boolean precedence.
    """
    return source.replace('|', ' or ').replace('&', ' and ')


def _replace_locals(source, local_symbol='@'):
    """Replace local variables with a syntacticall valid name."""
    return source.replace(local_symbol, _LOCAL_TAG)


def _preparse(source):
    """Compose assignment and boolean replacement."""
    return _replace_booleans(_rewrite_assign(source))


def _is_type(t):
    """Factory for a type checking function of type ``t`` or tuple of types."""
    return lambda x: isinstance(x.value, t)


_is_list = _is_type(list)
_is_str = _is_type(string_types)


# partition all AST nodes
_all_nodes = frozenset(filter(lambda x: isinstance(x, type) and
                              issubclass(x, ast.AST),
                              (getattr(ast, node) for node in dir(ast))))


def _filter_nodes(superclass, all_nodes=_all_nodes):
    """Filter out AST nodes that are subclasses of ``superclass``."""
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


_unsupported_expr_nodes = frozenset(['Yield', 'GeneratorExp', 'IfExp',
                                     'DictComp', 'SetComp', 'Repr', 'Lambda',
                                     'Set', 'AST', 'Is', 'IsNot'])

# these nodes are low priority or won't ever be supported (e.g., AST)
_unsupported_nodes = ((_stmt_nodes | _mod_nodes | _handler_nodes |
                       _arguments_nodes | _keyword_nodes | _alias_nodes |
                       _expr_context_nodes | _unsupported_expr_nodes) -
                      _hacked_nodes)

# we're adding a different assignment in some cases to be equality comparison
# and we don't want `stmt` and friends in their so get only the class whose
# names are capitalized
_base_supported_nodes = (_all_node_names - _unsupported_nodes) | _hacked_nodes
_msg = 'cannot both support and not support {0}'.format(_unsupported_nodes &
                                                        _base_supported_nodes)
assert not _unsupported_nodes & _base_supported_nodes, _msg


def _node_not_implemented(node_name, cls):
    """Return a function that raises a NotImplementedError with a passed node
    name.
    """

    def f(self, *args, **kwargs):
        raise NotImplementedError("{0!r} nodes are not "
                                  "implemented".format(node_name))
    return f


def disallow(nodes):
    """Decorator to disallow certain nodes from parsing. Raises a
    NotImplementedError instead.

    Returns
    -------
    disallowed : callable
    """
    def disallowed(cls):
        cls.unsupported_nodes = ()
        for node in nodes:
            new_method = _node_not_implemented(node, cls)
            name = 'visit_{0}'.format(node)
            cls.unsupported_nodes += (name,)
            setattr(cls, name, new_method)
        return cls
    return disallowed


def _op_maker(op_class, op_symbol):
    """Return a function to create an op class with its symbol already passed.

    Returns
    -------
    f : callable
    """

    def f(self, node, *args, **kwargs):
        """Return a partial function with an Op subclass with an operator
        already passed.

        Returns
        -------
        f : callable
        """
        return partial(op_class, op_symbol, *args, **kwargs)
    return f


_op_classes = {'binary': BinOp, 'unary': UnaryOp}


def add_ops(op_classes):
    """Decorator to add default implementation of ops."""
    def f(cls):
        for op_attr_name, op_class in compat.iteritems(op_classes):
            ops = getattr(cls, '{0}_ops'.format(op_attr_name))
            ops_map = getattr(cls, '{0}_op_nodes_map'.format(op_attr_name))
            for op in ops:
                op_node = ops_map[op]
                if op_node is not None:
                    made_op = _op_maker(op_class, op)
                    setattr(cls, 'visit_{0}'.format(op_node), made_op)
        return cls
    return f


@disallow(_unsupported_nodes)
@add_ops(_op_classes)
class BaseExprVisitor(ast.NodeVisitor):

    """Custom ast walker. Parsers of other engines should subclass this class
    if necessary.

    Parameters
    ----------
    env : Scope
    engine : str
    parser : str
    preparser : callable
    """
    const_type = Constant
    term_type = Term

    binary_ops = _cmp_ops_syms + _bool_ops_syms + _arith_ops_syms
    binary_op_nodes = ('Gt', 'Lt', 'GtE', 'LtE', 'Eq', 'NotEq', 'In', 'NotIn',
                       'BitAnd', 'BitOr', 'And', 'Or', 'Add', 'Sub', 'Mult',
                       None, 'Pow', 'FloorDiv', 'Mod')
    binary_op_nodes_map = dict(zip(binary_ops, binary_op_nodes))

    unary_ops = _unary_ops_syms
    unary_op_nodes = 'UAdd', 'USub', 'Invert', 'Not'
    unary_op_nodes_map = dict(zip(unary_ops, unary_op_nodes))

    rewrite_map = {
        ast.Eq: ast.In,
        ast.NotEq: ast.NotIn,
        ast.In: ast.In,
        ast.NotIn: ast.NotIn
    }

    def __init__(self, env, engine, parser, preparser=_preparse):
        self.env = env
        self.engine = engine
        self.parser = parser
        self.preparser = preparser
        self.assigner = None

    def visit(self, node, **kwargs):
        if isinstance(node, string_types):
            clean = self.preparser(node)
            node = ast.fix_missing_locations(ast.parse(clean))
        elif not isinstance(node, ast.AST):
            raise TypeError("Cannot visit objects of type {0!r}"
                            "".format(node.__class__.__name__))

        method = 'visit_' + node.__class__.__name__
        visitor = getattr(self, method)
        return visitor(node, **kwargs)

    def visit_Module(self, node, **kwargs):
        if len(node.body) != 1:
            raise SyntaxError('only a single expression is allowed')
        expr = node.body[0]
        return self.visit(expr, **kwargs)

    def visit_Expr(self, node, **kwargs):
        return self.visit(node.value, **kwargs)

    def _rewrite_membership_op(self, node, left, right):
        # the kind of the operator (is actually an instance)
        op_instance = node.op
        op_type = type(op_instance)

        # must be two terms and the comparison operator must be ==/!=/in/not in
        if is_term(left) and is_term(right) and op_type in self.rewrite_map:

            left_list, right_list = map(_is_list, (left, right))
            left_str, right_str = map(_is_str, (left, right))

            # if there are any strings or lists in the expression
            if left_list or right_list or left_str or right_str:
                op_instance = self.rewrite_map[op_type]()

            # pop the string variable out of locals and replace it with a list
            # of one string, kind of a hack
            if right_str:
                self.env.remove_tmp(right.name)
                name = self.env.add_tmp([right.value])
                right = self.term_type(name, self.env)

            if left_str:
                self.env.remove_tmp(left.name)
                name = self.env.add_tmp([left.value])
                left = self.term_type(name, self.env)

        op = self.visit(op_instance)
        return op, op_instance, left, right

    def _possibly_transform_eq_ne(self, node, left=None, right=None):
        if left is None:
            left = self.visit(node.left, side='left')
        if right is None:
            right = self.visit(node.right, side='right')
        op, op_class, left, right = self._rewrite_membership_op(node, left,
                                                                right)
        return op, op_class, left, right

    def _possibly_eval(self, binop, eval_in_python):
        # eval `in` and `not in` (for now) in "partial" python space
        # things that can be evaluated in "eval" space will be turned into
        # temporary variables. for example,
        # [1,2] in a + 2 * b
        # in that case a + 2 * b will be evaluated using numexpr, and the "in"
        # call will be evaluated using isin (in python space)
        return binop.evaluate(self.env, self.engine, self.parser,
                              self.term_type, eval_in_python)

    def _possibly_evaluate_binop(self, op, op_class, lhs, rhs,
                                 eval_in_python=('in', 'not in'),
                                 maybe_eval_in_python=('==', '!=', '<', '>',
                                                       '<=', '>=')):
        res = op(lhs, rhs)

        if self.engine != 'pytables':
            if (res.op in _cmp_ops_syms
                    and getattr(lhs, 'is_datetime', False)
                    or getattr(rhs, 'is_datetime', False)):
                # all date ops must be done in python bc numexpr doesn't work
                # well with NaT
                return self._possibly_eval(res, self.binary_ops)

        if res.op in eval_in_python:
            # "in"/"not in" ops are always evaluated in python
            return self._possibly_eval(res, eval_in_python)
        elif self.engine != 'pytables':
            if (getattr(lhs, 'return_type', None) == object
                    or getattr(rhs, 'return_type', None) == object):
                # evaluate "==" and "!=" in python if either of our operands
                # has an object return type
                return self._possibly_eval(res, eval_in_python +
                                           maybe_eval_in_python)
        return res

    def visit_BinOp(self, node, **kwargs):
        op, op_class, left, right = self._possibly_transform_eq_ne(node)
        return self._possibly_evaluate_binop(op, op_class, left, right)

    def visit_Div(self, node, **kwargs):
        return lambda lhs, rhs: Div(lhs, rhs,
                                    truediv=self.env.locals['truediv'])

    def visit_UnaryOp(self, node, **kwargs):
        op = self.visit(node.op)
        operand = self.visit(node.operand)
        return op(operand)

    def visit_Name(self, node, **kwargs):
        return self.term_type(node.id, self.env, **kwargs)

    def visit_Num(self, node, **kwargs):
        return self.const_type(node.n, self.env)

    def visit_Str(self, node, **kwargs):
        name = self.env.add_tmp(node.s)
        return self.term_type(name, self.env)

    def visit_List(self, node, **kwargs):
        name = self.env.add_tmp([self.visit(e).value for e in node.elts])
        return self.term_type(name, self.env)

    visit_Tuple = visit_List

    def visit_Index(self, node, **kwargs):
        """ df.index[4] """
        return self.visit(node.value)

    def visit_Subscript(self, node, **kwargs):
        value = self.visit(node.value)
        slobj = self.visit(node.slice)
        result = pd.eval(slobj, local_dict=self.env, engine=self.engine,
                         parser=self.parser)
        try:
            # a Term instance
            v = value.value[result]
        except AttributeError:
            # an Op instance
            lhs = pd.eval(value, local_dict=self.env, engine=self.engine,
                          parser=self.parser)
            v = lhs[result]
        name = self.env.add_tmp(v)
        return self.term_type(name, env=self.env)

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
        """
        support a single assignment node, like

        c = a + b

        set the assigner at the top level, must be a Name node which
        might or might not exist in the resolvers

        """

        if len(node.targets) != 1:
            raise SyntaxError('can only assign a single expression')
        if not isinstance(node.targets[0], ast.Name):
            raise SyntaxError('left hand side of an assignment must be a '
                              'single name')
        if self.env.target is None:
            raise ValueError('cannot assign without a target object')

        try:
            assigner = self.visit(node.targets[0], **kwargs)
        except UndefinedVariableError:
            assigner = node.targets[0].id

        self.assigner = getattr(assigner, 'name', assigner)
        if self.assigner is None:
            raise SyntaxError('left hand side of an assignment must be a '
                              'single resolvable name')

        return self.visit(node.value, **kwargs)

    def visit_Attribute(self, node, **kwargs):
        attr = node.attr
        value = node.value

        ctx = node.ctx
        if isinstance(ctx, ast.Load):
            # resolve the value
            resolved = self.visit(value).value
            try:
                v = getattr(resolved, attr)
                name = self.env.add_tmp(v)
                return self.term_type(name, self.env)
            except AttributeError:
                # something like datetime.datetime where scope is overridden
                if isinstance(value, ast.Name) and value.id == attr:
                    return resolved

        raise ValueError("Invalid Attribute context {0}".format(ctx.__name__))

    def visit_Call(self, node, side=None, **kwargs):

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
                raise ValueError("keyword error in function call "
                                 "'{0}'".format(node.func.id))
            keywords[key.arg] = self.visit(key.value).value
        if node.kwargs is not None:
            keywords.update(self.visit(node.kwargs).value)

        return self.const_type(res(*args, **keywords), self.env)

    def translate_In(self, op):
        return op

    def visit_Compare(self, node, **kwargs):
        ops = node.ops
        comps = node.comparators

        # base case: we have something like a CMP b
        if len(comps) == 1:
            op = self.translate_In(ops[0])
            binop = ast.BinOp(op=op, left=node.left, right=comps[0])
            return self.visit(binop)

        # recursive case: we have a chained comparison, a CMP b CMP c, etc.
        left = node.left
        values = []
        for op, comp in zip(ops, comps):
            new_node = self.visit(ast.Compare(comparators=[comp], left=left,
                                              ops=[self.translate_In(op)]))
            left = comp
            values.append(new_node)
        return self.visit(ast.BoolOp(op=ast.And(), values=values))

    def _try_visit_binop(self, bop):
        if isinstance(bop, (Op, Term)):
            return bop
        return self.visit(bop)

    def visit_BoolOp(self, node, **kwargs):
        def visitor(x, y):
            lhs = self._try_visit_binop(x)
            rhs = self._try_visit_binop(y)

            op, op_class, lhs, rhs = self._possibly_transform_eq_ne(node, lhs,
                                                                    rhs)
            return self._possibly_evaluate_binop(op, node.op, lhs, rhs)

        operands = node.values
        return reduce(visitor, operands)


_python_not_supported = frozenset(['Dict', 'Call', 'BoolOp', 'In', 'NotIn'])
_numexpr_supported_calls = frozenset(_reductions + _mathops)


@disallow((_unsupported_nodes | _python_not_supported) -
          (_boolop_nodes | frozenset(['BoolOp', 'Attribute', 'In', 'NotIn',
                                      'Tuple'])))
class PandasExprVisitor(BaseExprVisitor):

    def __init__(self, env, engine, parser,
                 preparser=lambda x: _replace_locals(_replace_booleans(x))):
        super(PandasExprVisitor, self).__init__(env, engine, parser, preparser)


@disallow(_unsupported_nodes | _python_not_supported | frozenset(['Not']))
class PythonExprVisitor(BaseExprVisitor):

    def __init__(self, env, engine, parser, preparser=lambda x: x):
        super(PythonExprVisitor, self).__init__(env, engine, parser,
                                                preparser=preparser)


class Expr(StringMixin):

    """Object encapsulating an expression.

    Parameters
    ----------
    expr : str
    engine : str, optional, default 'numexpr'
    parser : str, optional, default 'pandas'
    env : Scope, optional, default None
    truediv : bool, optional, default True
    level : int, optional, default 2
    """

    def __init__(self, expr, engine='numexpr', parser='pandas', env=None,
                 truediv=True, level=2):
        self.expr = expr
        self.env = _ensure_scope(level=level, local_dict=env)
        self.engine = engine
        self.parser = parser
        self._visitor = _parsers[parser](self.env, self.engine, self.parser)
        self.terms = self.parse()
        self.truediv = truediv

    @property
    def assigner(self):
        return getattr(self._visitor, 'assigner', None)

    def __call__(self):
        self.env.locals['truediv'] = self.truediv
        return self.terms(self.env)

    def __unicode__(self):
        return com.pprint_thing(self.terms)

    def __len__(self):
        return len(self.expr)

    def parse(self):
        """Parse an expression"""
        return self._visitor.visit(self.expr)

    def align(self):
        """align a set of Terms"""
        return self.terms.align(self.env)

    @property
    def names(self):
        """Get the names in an expression"""
        if is_term(self.terms):
            return frozenset([self.terms.name])
        return frozenset(term.name for term in com.flatten(self.terms))

    def check_name_clashes(self):
        env = self.env
        names = self.names
        res_keys = frozenset(env.resolver_dict.keys()) & names
        lcl_keys = frozenset(env.locals.keys()) & names
        gbl_keys = frozenset(env.globals.keys()) & names
        _check_disjoint_resolver_names(res_keys, lcl_keys, gbl_keys)

    def add_resolvers_to_locals(self):
        """Add the extra scope (resolvers) to local scope

        Notes
        -----
        This should be done after parsing and pre-evaluation, otherwise
        unnecessary name clashes will occur.
        """
        self.env.locals.update(self.env.resolver_dict)


def isexpr(s, check_names=True):
    """Strict checking for a valid expression."""
    try:
        Expr(s, env=_ensure_scope() if check_names else None)
    except SyntaxError:
        return False
    except NameError:
        return not check_names
    return True


_parsers = {'python': PythonExprVisitor, 'pandas': PandasExprVisitor}
