""" manage PyTables query interface via Expressions """

import ast
import time
import warnings
from functools import partial
from datetime import datetime

import numpy as np

import pandas.core.common as com
import pandas.lib as lib
from pandas.computation import expr, ops
from pandas.computation.ops import is_term, Constant
from pandas.computation.expr import BaseExprVisitor
from pandas import Index
from pandas.core.common import is_list_like


def _ensure_decoded(s):
    """ if we have bytes, decode them to unicode """
    if isinstance(s, (np.bytes_, bytes)):
        s = s.decode('UTF-8')
    return s


class Scope(expr.Scope):
    __slots__ = 'globals', 'locals', 'queryables'

    def __init__(self, gbls=None, lcls=None, queryables=None, frame_level=1):
        super(
            Scope,
            self).__init__(gbls=gbls,
                           lcls=lcls,
                           frame_level=frame_level)
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
        return self.env.locals.get(self.name,
                                   self.env.globals.get(self.name, self.name))


class BinOp(ops.BinOp):

    _max_selectors = 31

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
                if (isinstance(left, ConditionBinOp) and
                    isinstance(right, ConditionBinOp)):
                    k = JointConditionBinOp
                elif isinstance(left, k):
                    return left
                elif isinstance(right, k):
                    return right

            elif isinstance(left, FilterBinOp):
                if (isinstance(left, FilterBinOp) and
                    isinstance(right, FilterBinOp)):
                    k = JointFilterBinOp
                elif isinstance(left, k):
                    return left
                elif isinstance(right, k):
                    return right

            return k(self.op, left, right, queryables=self.queryables,
                     encoding=self.encoding).evaluate()

        left, right = self.lhs, self.rhs

        if is_term(left) and is_term(right):
            res = pr(left.value, right.value)
        elif not is_term(left) and is_term(right):
            res = pr(left.prune(klass), right.value)
        elif is_term(left) and not is_term(right):
            res = pr(left.value, right.prune(klass))
        elif not (is_term(left) or is_term(right)):
            res = pr(left.prune(klass), right.prune(klass))

        return res

    def conform(self, rhs):
        """ inplace conform rhs """
        if not is_list_like(rhs):
            rhs = [rhs]
        if hasattr(self.rhs, 'ravel'):
            rhs = rhs.ravel()
        return rhs

    @property
    def is_valid(self):
        """ return True if this is a valid field """
        return self.lhs in self.queryables

    @property
    def is_in_table(self):
        """ return True if this is a valid column name for generation (e.g. an
        actual column in the table) """
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
        """ convert the expression that is in the term to something that is
        accepted by pytables """

        def stringify(value):
            value = str(value)
            if self.encoding is not None:
                value = value.encode(self.encoding)
            return value

        kind = _ensure_decoded(self.kind)
        if kind == u'datetime64' or kind == u'datetime':

            if isinstance(v, (int, float)):
                v = stringify(v)
            v = _ensure_decoded(v)
            v = lib.Timestamp(v)
            if v.tz is not None:
                v = v.tz_convert('UTC')
            return TermValue(v, v.value, kind)
        elif isinstance(v, datetime) or hasattr(v, 'timetuple') or kind == u'date':
            v = time.mktime(v.timetuple())
            return TermValue(v, lib.Timestamp(v), kind)
        elif kind == u'integer':
            v = int(float(v))
            return TermValue(v, v, kind)
        elif kind == u'float':
            v = float(v)
            return TermValue(v, v, kind)
        elif kind == u'bool':
            if isinstance(v, basestring):
                v = not v.strip().lower() in [u'false', u'f', u'no', u'n',
                                              u'none', u'0', u'[]', u'{}', u'']
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
        return com.pprint_thing("[Filter : [{0}] -> "
                                "[{1}]".format(self.filter[0], self.filter[1]))

    def invert(self):
        """ invert the filter """
        if self.filter is not None:
            f = list(self.filter)
            f[1] = self.generate_filter_op(invert=True)
            self.filter = tuple(f)
        return self

    def format(self):
        """ return the actual filter format """
        return [self.filter]

    def evaluate(self):

        if not isinstance(self.lhs, basestring):
            return self

        if not self.is_valid:
            raise ValueError("query term is not valid [%s]" % self)

        rhs = self.conform(self.rhs)
        values = [TermValue(v, v, self.kind) for v in rhs]

        if self.is_in_table:

            # if too many values to create the expression, use a filter instead
            if self.op in ['==', '!='] and len(values) > self._max_selectors:

                filter_op = self.generate_filter_op()
                self.filter = (
                    self.lhs,
                    filter_op,
                    Index([v.value for v in values]))

                return self
            return None

        # equality conditions
        if self.op in ['==', '!=']:

            filter_op = self.generate_filter_op()
            self.filter = (
                self.lhs,
                filter_op,
                Index([v.value for v in values]))

        else:
            raise TypeError(
                "passing a filterable condition to a non-table indexer [%s]" %
                self)

        return self

    def generate_filter_op(self, invert=False):
        if (self.op == '!=' and not invert) or (self.op == '==' and invert):
            return lambda axis, vals: ~axis.isin(vals)
        else:
            return lambda axis, vals: axis.isin(vals)


class JointFilterBinOp(FilterBinOp):

    def format(self):
        raise NotImplementedError("unable to collapse Joint Filters")

    def evaluate(self):
        return self


class ConditionBinOp(BinOp):

    def __unicode__(self):
        return com.pprint_thing("[Condition : [{0}]]".format(self.condition))

    def invert(self):
        """ invert the condition """
        #if self.condition is not None:
        #    self.condition = "~(%s)" % self.condition
        #return self
        raise NotImplementedError("cannot use an invert condition when passing to numexpr")

    def format(self):
        """ return the actual ne format """
        return self.condition

    def evaluate(self):

        if not isinstance(self.lhs, basestring):
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
        self.condition = "(%s %s %s)" % (
            self.lhs.condition,
            self.op,
            self.rhs.condition)
        return self


class UnaryOp(ops.UnaryOp):

    def prune(self, klass):

        if self.op != '~':
            raise NotImplementedError("UnaryOp only support invert type ops")

        operand = self.operand
        operand = operand.prune(klass)

        if operand is not None:
            if issubclass(klass,ConditionBinOp):
                if operand.condition is not None:
                    return operand.invert()
            elif issubclass(klass,FilterBinOp):
                if operand.filter is not None:
                    return operand.invert()

        return None



_op_classes = {'unary': UnaryOp}


class ExprVisitor(BaseExprVisitor):
    def __init__(self, env, **kwargs):
        super(ExprVisitor, self).__init__(env)
        for bin_op in self.binary_ops:
            setattr(self, 'visit_{0}'.format(self.binary_op_nodes_map[bin_op]),
                    lambda node, bin_op=bin_op: partial(BinOp, bin_op,
                                                        **kwargs))

    def visit_Name(self, node, side=None, **kwargs):
        return Term(node.id, self.env, side=side, **kwargs)

    def visit_UnaryOp(self, node, **kwargs):
        if isinstance(node.op, (ast.Not, ast.Invert)):
            return UnaryOp('~', self.visit(node.operand))
        elif isinstance(node.op, ast.USub):
            return Constant(-self.visit(node.operand).value, self.env)
        elif isinstance(node.op, ast.UAdd):
            raise NotImplementedError('Unary addition not supported')

    def visit_USub(self, node, **kwargs):
        return Constant(-self.visit(node.operand).value, self.env)

    def visit_Index(self, node, **kwargs):
        return self.visit(node.value).value

class Expr(expr.Expr):

    """ hold a pytables like expression, comprised of possibly multiple 'terms'

    Parameters
    ----------
    where : string term expression, Expr, or list-like of Exprs
    queryables : a kinds map (dict of column name -> kind), or None i column is non-indexable
    encoding : an encoding that will encode the query terms

    Returns
    -------
    an Expr object

    Examples
    --------

    'index>=date'
    "columns=['A', 'D']"
    'columns=A'
    'columns==A'
    "~(columns=['A','B'])"
    'index>df.index[3] & string="bar"'
    '(index>df.index[3] & index<=df.index[6]) | string="bar"'
    "ts>=Timestamp('2012-02-01')"
    "major_axis>=20130101"
    """

    def __init__(self, where, op=None, value=None, queryables=None,
                 encoding=None, scope_level=None):

        # try to be back compat
        where = self.parse_back_compat(where, op, value)

        self.encoding = encoding
        self.condition = None
        self.filter = None
        self.terms = None
        self._visitor = None

        # capture the environement if needed
        lcls = dict()
        if isinstance(where, Expr):

            lcls.update(where.env.locals)
            where = str(where)

        elif isinstance(where, (list, tuple)):

            for w in where:
                if isinstance(w, Expr):
                    lcls.update(w.env.locals)
                else:
                    w = self.parse_back_compat(w)

            where = ' & ' .join(["(%s)" % w for w in where])

        self.expr = where
        self.env = Scope(lcls=lcls)
        self.env.update(scope_level)

        if queryables is not None:
            self.env.queryables.update(queryables)
            self._visitor = ExprVisitor(self.env, queryables=queryables,
                                        encoding=encoding)
            self.terms = self.parse()

    def parse_back_compat(self, w, op=None, value=None):
        """ allow backward compatibility for passed arguments """

        if isinstance(w, dict):
            w, op, value = w.get('field'), w.get('op'), w.get('value')
            if not isinstance(w, basestring):
                raise TypeError(
                    "where must be passed as a string if op/value are passed")
            warnings.warn("passing a dict to Expr is deprecated, "
                          "pass the where as a single string",
                          DeprecationWarning)

        if op is not None:
            if not isinstance(w, basestring):
                raise TypeError(
                    "where must be passed as a string if op/value are passed")

            if isinstance(op, Expr):
                raise TypeError("invalid op passed, must be a string")
            w = "{0}{1}".format(w, op)
            if value is not None:
                if isinstance(value, Expr):
                    raise TypeError("invalid value passed, must be a string")
                w = "{0}{1}".format(w, value)

            warnings.warn("passing multiple values to Expr is deprecated, "
                          "pass the where as a single string",
                          DeprecationWarning)

        return w

    def __unicode__(self):
        if self.terms is not None:
            return unicode(self.terms)
        return self.expr

    def evaluate(self):
        """ create and return the numexpr condition and filter """

        try:
            self.condition = self.terms.prune(ConditionBinOp)
        except AttributeError:
            raise ValueError(
                "cannot process expression [{0}], [{1}] is not a valid condition".format(self.expr,self))
        try:
            self.filter = self.terms.prune(FilterBinOp)
        except AttributeError:
            raise ValueError(
                "cannot process expression [{0}], [{1}] is not a valid filter".format(self.expr,self))

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
