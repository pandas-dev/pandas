import sys
import re
import ast
from functools import partial

from pandas.computation import expr, ops
from pandas.computation.ops import is_term
from pandas.computation.expr import ExprParserError

class Scope(expr.Scope):
    __slots__ = 'globals', 'locals', 'queryables'

    def __init__(self, gbls=None, lcls=None, queryables=None, frame_level=1):
        super(Scope, self).__init__(gbls=gbls, lcls=lcls, frame_level=frame_level)
        self.queryables = queryables or dict()

class Term(ops.Term):

    def __init__(self, name, env, side=None):
        super(Term, self).__init__(name, env, side=side)

    def _resolve_name(self):

        # must be a queryable
        if self.side == 'left':
            if self.name not in self.env.queryables:
                raise NameError('name {0!r} is not defined'.format(self.name))
            return self.name

        # resolve the rhs (and allow to be None)
        return self.env.locals.get(self.name, self.env.globals.get(self.name,self.name))

def format_value(q, lhs, v):
    """ given a queryable, a lhs name and value, return a formatted value """
    return v

class BinOp(ops.BinOp):

    def __call__(self, q):
        left, right = self.lhs, self.rhs

        # base cases
        if is_term(left) and is_term(right):
            res = "(%s %s %s)" % (left.value,self.op,format_value(q, left.value, right.value))
        elif not is_term(left) and is_term(right):
            res = "(%s %s %s)" % (left(q),self.op,right.value)
        elif is_term(left) and not is_term(right):
            res = "(%s %s %s)" % (left.value,self.op,right(q))
        elif not (is_term(left) or is_term(right)):
            res = "(%s %s %s)" % (left(q),self.op,right(q))

        return res

class UnaryOp(ops.UnaryOp):
    def __call__(self, q):
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

    def __init__(self, env):
        for bin_op in self.bin_ops:
            setattr(self, 'visit_{0}'.format(self.bin_op_nodes_map[bin_op]),
                    lambda node, bin_op=bin_op: partial(BinOp, bin_op))

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

    _max_selectors = 31

    def __init__(self, expression, queryables=None, encoding=None):
        self.expr = self.pre_parse(expression)
        self.env = Scope(queryables=queryables,frame_level=2)
        self._visitor = ExprVisitor(self.env)
        self.terms = self.parse()
        self.encoding = encoding
        self.condition = None
        self.filter = None

    def pre_parse(self, expression):
        """ transform = to == """
        expression = re.sub("=+","==",expression)
        return expression

    def evaluate(self):
        """ create and return the numexpr condition and filter """
        import pdb; pdb.set_trace()
        terms = []
        filter = []

        self.terms(self.env)
        #for t in self.terms:

        terms = [t for t in self.terms if t.condition is not None]
        if len(terms):
            self.condition = "(%s)" % ' & '.join(
                [t.condition for t in terms])
            self.filter = []
            for t in self.terms:
                if t.filter is not None:
                    self.filter.append(t.filter)


    @property
    def is_valid(self):
        """ return True if this is a valid field """
        return self.field in self.q

    @property
    def is_in_table(self):
        """ return True if this is a valid column name for generation (e.g. an actual column in the table) """
        return self.q.get(self.field) is not None

    @property
    def kind(self):
        """ the kind of my field """
        return self.q.get(self.field)

    def generate(self, v):
        """ create and return the op string for this TermValue """
        val = v.tostring(self.encoding)
        return "(%s %s %s)" % (self.field, self.op, val)

        """ set the numexpr expression for this term """

        if not self.is_valid:
            raise ValueError("query term is not valid [%s]" % str(self))

        # convert values if we are in the table
        if self.is_in_table:
            values = [self.convert_value(v) for v in self.value]
        else:
            values = [TermValue(v, v, self.kind) for v in self.value]

        # equality conditions
        if self.op in ['==', '!=']:

            # our filter op expression
            if self.op == '!=':
                filter_op = lambda axis, vals: not axis.isin(vals)
            else:
                filter_op = lambda axis, vals: axis.isin(vals)

            if self.is_in_table:

                # too many values to create the expression?
                if len(values) <= self._max_selectors:
                    vs = [self.generate(v) for v in values]
                    self.condition = "(%s)" % ' | '.join(vs)

                # use a filter after reading
                else:
                    self.filter = (
                        self.field,
                        filter_op,
                        Index([v.value for v in values]))

            else:

                self.filter = (
                    self.field,
                    filter_op,
                    Index([v.value for v in values]))

        else:

            if self.is_in_table:

                self.condition = self.generate(values[0])

            else:

                raise TypeError(
                    "passing a filterable condition to a non-table indexer [%s]" %
                    str(self))

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


