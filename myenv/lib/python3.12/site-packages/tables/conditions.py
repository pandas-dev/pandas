"""Utility functions and classes for supporting query conditions.

Classes:

`CompileCondition`
    Container for a compiled condition.

Functions:

`compile_condition`
    Compile a condition and extract usable index conditions.
`call_on_recarr`
    Evaluate a function over a structured array.

"""

import re
import numexpr as ne

from .utilsextension import get_nested_field
from .utils import lazyattr


_no_matching_opcode = re.compile(r"[^a-z]([a-z]+)_([a-z]+)[^a-z]")
# E.g. "gt" and "bfc" from "couldn't find matching opcode for 'gt_bfc'".


def _unsupported_operation_error(exception):
    """Make the \"no matching opcode\" Numexpr `exception` more clear.

    A new exception of the same kind is returned.

    """

    message = exception.args[0]
    op, types = _no_matching_opcode.search(message).groups()
    newmessage = "unsupported operand types for *%s*: " % op
    newmessage += ', '.join(
        ne.necompiler.typecode_to_kind[t] for t in types[1:])
    return exception.__class__(newmessage)


def _check_indexable_cmp(getidxcmp):
    """Decorate `getidxcmp` to check the returned indexable comparison.

    This does some extra checking that Numexpr would perform later on
    the comparison if it was compiled within a complete condition.

    """

    def newfunc(exprnode, indexedcols):
        result = getidxcmp(exprnode, indexedcols)
        if result[0] is not None:
            try:
                ne.necompiler.typeCompileAst(
                    ne.necompiler.expressionToAST(exprnode))
            except NotImplementedError as nie:
                # Try to make this Numexpr error less cryptic.
                raise _unsupported_operation_error(nie)
        return result
    newfunc.__name__ = getidxcmp.__name__
    newfunc.__doc__ = getidxcmp.__doc__
    return newfunc


@_check_indexable_cmp
def _get_indexable_cmp(exprnode, indexedcols):
    """Get the indexable variable-constant comparison in `exprnode`.

    A tuple of (variable, operation, constant) is returned if
    `exprnode` is a variable-constant (or constant-variable)
    comparison, and the variable is in `indexedcols`.  A normal
    variable can also be used instead of a constant: a tuple with its
    name will appear instead of its value.

    Otherwise, the values in the tuple are ``None``.
    """

    not_indexable = (None, None, None)
    turncmp = {'lt': 'gt',
               'le': 'ge',
               'eq': 'eq',
               'ge': 'le',
               'gt': 'lt', }

    def get_cmp(var, const, op):
        var_value, const_value = var.value, const.value
        if (var.astType == 'variable' and var_value in indexedcols
           and const.astType in ['constant', 'variable']):
            if const.astType == 'variable':
                const_value = (const_value, )
            return (var_value, op, const_value)
        return None

    def is_indexed_boolean(node):
        return (node.astType == 'variable'
                and node.astKind == 'bool'
                and node.value in indexedcols)

    # Boolean variables are indexable by themselves.
    if is_indexed_boolean(exprnode):
        return (exprnode.value, 'eq', True)
    # And so are negations of boolean variables.
    if exprnode.astType == 'op' and exprnode.value == 'invert':
        child = exprnode.children[0]
        if is_indexed_boolean(child):
            return (child.value, 'eq', False)
        # A negation of an expression will be returned as ``~child``.
        # The indexability of the negated expression will be decided later on.
        if child.astKind == "bool":
            return (child, 'invert', None)

    # Check node type.  Only comparisons are indexable from now on.
    if exprnode.astType != 'op':
        return not_indexable
    cmpop = exprnode.value
    if cmpop not in turncmp:
        return not_indexable

    # Look for a variable-constant comparison in both directions.
    left, right = exprnode.children
    cmp_ = get_cmp(left, right, cmpop)
    if cmp_:
        return cmp_
    cmp_ = get_cmp(right, left, turncmp[cmpop])
    if cmp_:
        return cmp_

    return not_indexable


def _equiv_expr_node(x, y):
    """Returns whether two ExpressionNodes are equivalent.

    This is needed because '==' is overridden on ExpressionNode to
    return a new ExpressionNode.

    """
    if (not isinstance(x, ne.expressions.ExpressionNode)
            and not isinstance(y, ne.expressions.ExpressionNode)):
        return x == y
    elif (type(x) is not type(y)
          or not isinstance(x, ne.expressions.ExpressionNode)
          or not isinstance(y, ne.expressions.ExpressionNode)
          or x.value != y.value
          or x.astKind != y.astKind
          or len(x.children) != len(y.children)):
        return False
    for xchild, ychild in zip(x.children, y.children):
        if not _equiv_expr_node(xchild, ychild):
            return False
    return True


def _get_idx_expr_recurse(exprnode, indexedcols, idxexprs, strexpr):
    """Here lives the actual implementation of the get_idx_expr() wrapper.

    'idxexprs' is a list of expressions in the form ``(var, (ops),
    (limits))``. 'strexpr' is the indexable expression in string format.
    These parameters will be received empty (i.e. [], ['']) for the
    first time and populated during the different recursive calls.
    Finally, they are returned in the last level to the original
    wrapper.  If 'exprnode' is not indexable, it will return the tuple
    ([], ['']) so as to signal this.

    """

    not_indexable = ([], [''])
    op_conv = {
        'and': '&',
        'or': '|',
        'not': '~',
    }
    negcmp = {
        'lt': 'ge',
        'le': 'gt',
        'ge': 'lt',
        'gt': 'le',
    }

    def fix_invert(idxcmp, exprnode, indexedcols):
        invert = False
        # Loop until all leading negations have been dealt with
        while idxcmp[1] == "invert":
            invert ^= True
            # The information about the negated node is in first position
            exprnode = idxcmp[0]
            idxcmp = _get_indexable_cmp(exprnode, indexedcols)
        return idxcmp, exprnode, invert

    # Indexable variable-constant comparison.
    idxcmp = _get_indexable_cmp(exprnode, indexedcols)
    idxcmp, exprnode, invert = fix_invert(idxcmp, exprnode, indexedcols)
    if idxcmp[0]:
        if invert:
            var, op, value = idxcmp
            if op == 'eq' and value in [True, False]:
                # ``var`` must be a boolean index.  Flip its value.
                value ^= True
            else:
                op = negcmp[op]
            expr = (var, (op,), (value,))
            invert = False
        else:
            expr = (idxcmp[0], (idxcmp[1],), (idxcmp[2],))
        return [expr]

    # For now negations of complex expressions will be not supported as
    # forming part of an indexable condition.  This might be supported in
    # the future.
    if invert:
        return not_indexable

    # Only conjunctions and disjunctions of comparisons are considered
    # for the moment.
    if exprnode.astType != 'op' or exprnode.value not in ['and', 'or']:
        return not_indexable

    left, right = exprnode.children
    # Get the expression at left
    lcolvar, lop, llim = _get_indexable_cmp(left, indexedcols)
    # Get the expression at right
    rcolvar, rop, rlim = _get_indexable_cmp(right, indexedcols)

    # Use conjunction of indexable VC comparisons like
    # ``(a <[=] x) & (x <[=] b)`` or ``(a >[=] x) & (x >[=] b)``
    # as ``a <[=] x <[=] b``, for the moment.
    op = exprnode.value
    if (lcolvar is not None and rcolvar is not None
            and _equiv_expr_node(lcolvar, rcolvar) and op == 'and'):
        if lop in ['gt', 'ge'] and rop in ['lt', 'le']:  # l <= x <= r
            expr = (lcolvar, (lop, rop), (llim, rlim))
            return [expr]
        if lop in ['lt', 'le'] and rop in ['gt', 'ge']:  # l >= x >= r
            expr = (rcolvar, (rop, lop), (rlim, llim))
            return [expr]

    # Recursively get the expressions at the left and the right
    lexpr = _get_idx_expr_recurse(left, indexedcols, idxexprs, strexpr)
    rexpr = _get_idx_expr_recurse(right, indexedcols, idxexprs, strexpr)

    def add_expr(expr, idxexprs, strexpr):
        """Add a single expression to the list."""

        if isinstance(expr, list):
            # expr is a single expression
            idxexprs.append(expr[0])
            lenexprs = len(idxexprs)
            # Mutate the strexpr string
            if lenexprs == 1:
                strexpr[:] = ["e0"]
            else:
                strexpr[:] = [
                    "(%s %s e%d)" % (strexpr[0], op_conv[op], lenexprs - 1)]

    # Add expressions to the indexable list when they are and'ed, or
    # they are both indexable.
    if lexpr != not_indexable and (op == "and" or rexpr != not_indexable):
        add_expr(lexpr, idxexprs, strexpr)
        if rexpr != not_indexable:
            add_expr(rexpr, idxexprs, strexpr)
        return (idxexprs, strexpr)
    if rexpr != not_indexable and op == "and":
        add_expr(rexpr, idxexprs, strexpr)
        return (idxexprs, strexpr)

    # Can not use indexed column.
    return not_indexable


def _get_idx_expr(expr, indexedcols):
    """Extract an indexable expression out of `exprnode`.

    Looks for variable-constant comparisons in the expression node
    `exprnode` involving variables in `indexedcols`.

    It returns a tuple of (idxexprs, strexpr) where 'idxexprs' is a
    list of expressions in the form ``(var, (ops), (limits))`` and
    'strexpr' is the indexable expression in string format.

    Expressions such as ``0 < c1 <= 1`` do not work as expected.

    Right now only some of the *indexable comparisons* are considered:

    * ``a <[=] x``, ``a == x`` and ``a >[=] x``
    * ``(a <[=] x) & (y <[=] b)`` and ``(a == x) | (b == y)``
    * ``~(~c_bool)``, ``~~c_bool`` and ``~(~c_bool) & (c_extra != 2)``

    (where ``a``, ``b`` and ``c_bool`` are indexed columns, but
    ``c_extra`` is not)

    Particularly, the ``!=`` operator and negations of complex boolean
    expressions are *not considered* as valid candidates:

    * ``a != 1`` and  ``c_bool != False``
    * ``~((a > 0) & (c_bool))``

    """

    return _get_idx_expr_recurse(expr, indexedcols, [], [''])


class CompiledCondition:
    """Container for a compiled condition."""

    @lazyattr
    def index_variables(self):
        """The columns participating in the index expression."""

        idxexprs = self.index_expressions
        idxvars = []
        for expr in idxexprs:
            idxvar = expr[0]
            if idxvar not in idxvars:
                idxvars.append(idxvar)
        return frozenset(idxvars)

    def __init__(self, func, params, idxexprs, strexpr, **kwargs):
        self.function = func
        """The compiled function object corresponding to this condition."""
        self.parameters = params
        """A list of parameter names for this condition."""
        self.index_expressions = idxexprs
        """A list of expressions in the form ``(var, (ops), (limits))``."""
        self.string_expression = strexpr
        """The indexable expression in string format."""
        self.kwargs = kwargs
        """NumExpr kwargs (used to pass ex_uses_vml to numexpr)"""

    def __repr__(self):
        return ("idxexprs: %s\nstrexpr: %s\nidxvars: %s"
                % (self.index_expressions, self.string_expression,
                   self.index_variables))

    def with_replaced_vars(self, condvars):
        """Replace index limit variables with their values in-place.

        A new compiled condition is returned.  Values are taken from
        the `condvars` mapping and converted to Python scalars.
        """

        exprs = self.index_expressions
        exprs2 = []
        for expr in exprs:
            idxlims = expr[2]  # the limits are in third place
            limit_values = []
            for idxlim in idxlims:
                if isinstance(idxlim, tuple):  # variable
                    idxlim = condvars[idxlim[0]]  # look up value
                    idxlim = idxlim.tolist()  # convert back to Python
                limit_values.append(idxlim)
            # Add this replaced entry to the new exprs2
            var, ops, _ = expr
            exprs2.append((var, ops, tuple(limit_values)))
        # Create a new container for the converted values
        newcc = CompiledCondition(
            self.function, self.parameters, exprs2, self.string_expression,
            **self.kwargs)
        return newcc


def _get_variable_names(expression):
    """Return the list of variable names in the Numexpr `expression`."""

    names = []
    stack = [expression]
    while stack:
        node = stack.pop()
        if node.astType == 'variable':
            names.append(node.value)
        elif hasattr(node, 'children'):
            stack.extend(node.children)
    return list(set(names))  # remove repeated names


def compile_condition(condition, typemap, indexedcols):
    """Compile a condition and extract usable index conditions.

    Looks for variable-constant comparisons in the `condition` string
    involving the indexed columns whose variable names appear in
    `indexedcols`.  The part of `condition` having usable indexes is
    returned as a compiled condition in a `CompiledCondition` container.

    Expressions such as '0 < c1 <= 1' do not work as expected.  The
    Numexpr types of *all* variables must be given in the `typemap`
    mapping.  The ``function`` of the resulting `CompiledCondition`
    instance is a Numexpr function object, and the ``parameters`` list
    indicates the order of its parameters.

    """

    # Get the expression tree and extract index conditions.
    expr = ne.necompiler.stringToExpression(condition, typemap, {})
    if expr.astKind != 'bool':
        raise TypeError("condition ``%s`` does not have a boolean type"
                        % condition)
    idxexprs = _get_idx_expr(expr, indexedcols)
    # Post-process the answer
    if isinstance(idxexprs, list):
        # Simple expression
        strexpr = ['e0']
    else:
        # Complex expression
        idxexprs, strexpr = idxexprs
    # Get rid of the unneccessary list wrapper for strexpr
    strexpr = strexpr[0]

    # Get the variable names used in the condition.
    # At the same time, build its signature.
    varnames = _get_variable_names(expr)
    signature = [(var, typemap[var]) for var in varnames]
    try:
        # See the comments in `numexpr.evaluate()` for the
        # reasons of inserting copy operators for unaligned,
        # *unidimensional* arrays.
        func = ne.necompiler.NumExpr(expr, signature)
    except NotImplementedError as nie:
        # Try to make this Numexpr error less cryptic.
        raise _unsupported_operation_error(nie)

    _, ex_uses_vml = ne.necompiler.getExprNames(condition, {})
    kwargs = {'ex_uses_vml': ex_uses_vml}

    params = varnames
    # This is more comfortable to handle about than a tuple.
    return CompiledCondition(func, params, idxexprs, strexpr, **kwargs)


def call_on_recarr(func, params, recarr, param2arg=None, **kwargs):
    """Call `func` with `params` over `recarr`.

    The `param2arg` function, when specified, is used to get an argument
    given a parameter name; otherwise, the parameter itself is used as
    an argument.  When the argument is a `Column` object, the proper
    column from `recarr` is used as its value.

    """

    args = []
    for param in params:
        if param2arg:
            arg = param2arg(param)
        else:
            arg = param
        if hasattr(arg, 'pathname'):  # looks like a column
            arg = get_nested_field(recarr, arg.pathname)
        args.append(arg)
    return func(*args, **kwargs)
