import sys, inspect
import re
import ast
from functools import partial
import numpy as np
from datetime import datetime
import time

import pandas
import pandas.core.common as com
import pandas.lib as lib
from pandas.computation import expr, ops
from pandas.computation.ops import is_term, Value
from pandas.computation.expr import ExprParserError
from pandas import Index
from pandas.core.common import is_list_like

def _ensure_decoded(s):
    """ if we have bytes, decode them to unicde """
    if isinstance(s, np.bytes_):
        s = s.decode('UTF-8')
    return s

class Scope(expr.Scope):
    __slots__ = 'globals', 'locals', 'queryables'

    def __init__(self, gbls=None, lcls=None, queryables=None, frame_level=1):
        super(Scope, self).__init__(gbls=gbls, lcls=lcls, frame_level=frame_level)
        self.queryables = queryables or dict()

class Term(ops.Term):

    def __init__(self, name, env, side=None):
        super(Term, self).__init__(name, env, side=side)

    def _resolve_name(self):

        # must be a queryables
        if self.side == 'left':
            if self.name not in self.env.queryables:
                raise NameError('name {0!r} is not defined'.format(self.name))
            return self.name

        # resolve the rhs (and allow to be None)
        return self.env.locals.get(self.name, self.env.globals.get(self.name,self.name))

class BinOp(ops.BinOp):

    def __init__(self, op, lhs, rhs, queryables, encoding):
        super(BinOp, self).__init__(op, lhs, rhs)
        self.queryables = queryables
        self.encoding = encoding
        self.filter = None
        self.condition = None

    def prune(self, klass):

        def pr(left, right):
            """ create and return a new specilized BinOp from myself """

            if left is None:
                return right
            elif right is None:
                return left

            k = klass
            if isinstance(left, ConditionBinOp):
                if isinstance(left, ConditionBinOp) and isinstance(right, ConditionBinOp):
                    k = JointConditionBinOp
                elif isinstance(left, k):
                    return left
                elif isinstance(right, k):
                    return right

            elif isinstance(left, FilterBinOp):
                if isinstance(left, FilterBinOp) and isinstance(right, FilterBinOp):
                    k = JointFilterBinOp
                elif isinstance(left, k):
                    return left
                elif isinstance(right, k):
                    return right

            return k(self.op, left, right, queryables=self.queryables, encoding=self.encoding).evaluate()

        left, right = self.lhs, self.rhs

        if is_term(left) and is_term(right):
            res = pr(left.value,right.value)
        elif not is_term(left) and is_term(right):
            res = pr(left.prune(klass),right.value)
        elif is_term(left) and not is_term(right):
            res = pr(left.value,right.prune(klass))
        elif not (is_term(left) or is_term(right)):
            res = pr(left.prune(klass),right.prune(klass))

        return res

    def conform(self, rhs):
        """ inplace conform rhs """
        if not is_list_like(rhs):
            rhs = [ rhs ]
        if hasattr(self.rhs,'ravel'):
            rhs = rhs.ravel()
        return rhs

    @property
    def is_valid(self):
        """ return True if this is a valid field """
        return self.lhs in self.queryables

    @property
    def is_in_table(self):
        """ return True if this is a valid column name for generation (e.g. an actual column in the table) """
        return self.queryables.get(self.lhs) is not None

    @property
    def kind(self):
        """ the kind of my field """
        return self.queryables.get(self.lhs)

    def generate(self, v):
        """ create and return the op string for this TermValue """
        val = v.tostring(self.encoding)
        return "(%s %s %s)" % (self.lhs, self.op, val)

    def convert_value(self, v):
        """ convert the expression that is in the term to something that is accepted by pytables """

        def stringify(value):
            value = str(value)
            if self.encoding is not None:
                value = value.encode(self.encoding)
            return value

        kind = _ensure_decoded(self.kind)
        if kind == u'datetime64' or kind == u'datetime':
            v = lib.Timestamp(v)
            if v.tz is not None:
                v = v.tz_convert('UTC')
            return TermValue(v, v.value, kind)
        elif isinstance(v, datetime) or hasattr(v, 'timetuple') or kind == u'date':
            v = time.mktime(v.timetuple())
            return TermValue(v, Timestamp(v), kind)
        elif kind == u'integer':
            v = int(float(v))
            return TermValue(v, v, kind)
        elif kind == u'float':
            v = float(v)
            return TermValue(v, v, kind)
        elif kind == u'bool':
            if isinstance(v, basestring):
                v = not v.strip().lower() in [
                    u'false', u'f', u'no', u'n', u'none', u'0', u'[]', u'{}', u'']
            else:
                v = bool(v)
            return TermValue(v, v, kind)
        elif not isinstance(v, basestring):
            v = stringify(v)
            return TermValue(v, stringify(v), u'string')

        # string quoting
        return TermValue(v, stringify(v), u'string')

class FilterBinOp(BinOp):

    def __unicode__(self):
        return com.pprint_thing("[Filter : [{0}] -> [{1}]".format(self.filter[0],self.filter[1]))

    def format(self):
        """ return the actual filter format """
        return [ self.filter ]

    def evaluate(self):

        if not isinstance(self.lhs,basestring):
            return self

        if not self.is_valid:
            raise ValueError("query term is not valid [%s]" % self)

        if self.is_in_table:
            return None

        rhs = self.conform(self.rhs)
        values = [TermValue(v, v, self.kind) for v in rhs]

        # equality conditions
        if self.op in ['==', '!=']:

            # our filter op expression
            if self.op == '!=':
                filter_op = lambda axis, vals: not axis.isin(vals)
            else:
                filter_op = lambda axis, vals: axis.isin(vals)

            self.filter = (
                self.lhs,
                filter_op,
                Index([v.value for v in values]))

        else:

            raise TypeError(
                "passing a filterable condition to a non-table indexer [%s]" %
                self)

        return self

class JointFilterBinOp(FilterBinOp):

    def format(self):
        raise NotImplementedError("unable to collapse Joint Filters")

    def evaluate(self):
        return self

class ConditionBinOp(BinOp):

    _max_selectors = 31

    def __unicode__(self):
        return com.pprint_thing("[Condition : [{0}]]".format(self.condition))

    def format(self):
        """ return the actual ne format """
        return self.condition

    def evaluate(self):

        if not isinstance(self.lhs,basestring):
            return self

        if not self.is_valid:
            raise ValueError("query term is not valid [%s]" % self)

        # convert values if we are in the table
        if not self.is_in_table:
            return None

        rhs = self.conform(self.rhs)
        values = [self.convert_value(v) for v in rhs]

        # equality conditions
        if self.op in ['==', '!=']:

            # our filter op expression
            if self.op == '!=':
                filter_op = lambda axis, vals: not axis.isin(vals)
            else:
                filter_op = lambda axis, vals: axis.isin(vals)

            # too many values to create the expression?
            if len(values) <= self._max_selectors:
                vs = [self.generate(v) for v in values]
                self.condition = "(%s)" % ' | '.join(vs)

            # use a filter after reading
            else:
                return None

        else:

            self.condition = self.generate(values[0])

        return self

class JointConditionBinOp(ConditionBinOp):

    def evaluate(self):
        self.condition = "(%s %s %s)" % (self.lhs.condition,self.op,self.rhs.condition)
        return self

class UnaryOp(ops.UnaryOp):

    def apply(self, func):
        operand = self.operand
        v = operand.value if is_term(operand) else operand
        return "%s (%s)" % (operand,v)

class ExprVisitor(expr.ExprVisitor):

    bin_ops = '>', '<', '>=', '<=', '==', '!=', '&', '|'
    bin_op_nodes = ('Gt', 'Lt', 'GtE', 'LtE', 'Eq', 'NotEq', 'BitAnd', 'BitOr')
    bin_op_nodes_map = dict(zip(bin_ops, bin_op_nodes))

    unary_ops =  ['~']
    unary_op_nodes = 'Invert'
    unary_op_nodes_map = dict(zip(unary_ops, unary_op_nodes))

    def __init__(self, env, **kwargs):
        for bin_op in self.bin_ops:
            setattr(self, 'visit_{0}'.format(self.bin_op_nodes_map[bin_op]),
                    lambda node, bin_op=bin_op: partial(BinOp, bin_op, **kwargs))

        for unary_op in self.unary_ops:
            setattr(self,
                    'visit_{0}'.format(self.unary_op_nodes_map[unary_op]),
                    lambda node, unary_op=unary_op: partial(UnaryOp, unary_op))
        self.env = env

    def visit_Module(self, node, **kwargs):
        if len(node.body) != 1:
            raise ExprParserError('only a single expression is allowed')

        body = node.body[0]
        return self.visit(body)

    def visit_Attribute(self, node, **kwargs):
        attr = node.attr
        value = node.value

        # resolve the value
        return getattr(self.visit(value).value,attr)

    def visit_Call(self, node, **kwargs):
        if not isinstance(node.func, ast.Name):
            raise TypeError("Only named functions are supported")

        res = self.visit(node.func)
        if res is None:
            raise ValueError("Invalid function call {0}".format(node.func.id))
        return res

    def visit_Compare(self, node, **kwargs):
        ops = node.ops
        comps = node.comparators
        for op, comp in zip(ops, comps):
            node = self.visit(op)(self.visit(node.left,side='left'), self.visit(comp,side='right'))
        return node

    def visit_Name(self, node, side=None, **kwargs):
        return Term(node.id, self.env, side=side)

class Expr(expr.Expr):

    """ hold a pytables like expression, comprised of possibly multiple 'terms'

    Parameters
    ----------
    field : dict, string term expression, or the field to operate (must be a valid index/column type of DataFrame/Panel)
    queryables : a kinds map (dict of column name -> kind), or None i column is non-indexable
    encoding : an encoding that will encode the query terms

    Returns
    -------
    an Expr object

    Examples
    --------
    """

    def __init__(self, expression, queryables=None, encoding=None, lcls=None):
        if isinstance(expression, Expr):
            expression = str(expression)
        self.expr = expression
        self.condition = None
        self.filter = None
        self.terms = None
        self._visitor = None

        # add current locals scope
        frame = inspect.currentframe()
        try:
            if lcls is None:
                lcls = dict()
            lcls.update(frame.f_back.f_locals)
            self.env = Scope(lcls = lcls)
        finally:
            del frame

        if queryables is not None:

            # if using the old format, this will raise
            if not isinstance(queryables, dict):
                raise TypeError("Expr must be called with a single-string expression")

            self.env.queryables.update(queryables)
            self._visitor = ExprVisitor(self.env, queryables=queryables, encoding=encoding)
            self.terms = self.parse()

    def __unicode__(self):
        if self.terms is not None:
            return unicode(self.terms)
        return self.expr

    def evaluate(self):
        """ create and return the numexpr condition and filter """

        self.condition = self.terms.prune(ConditionBinOp)
        self.filter = self.terms.prune(FilterBinOp)

        return self.condition, self.filter

class TermValue(object):

    """ hold a term value the we use to construct a condition/filter """

    def __init__(self, value, converted, kind):
        self.value = value
        self.converted = converted
        self.kind = kind

    def tostring(self, encoding):
        """ quote the string if not encoded
            else encode and return """
        if self.kind == u'string':
            if encoding is not None:
                return self.converted
            return '"%s"' % self.converted
        return self.converted


