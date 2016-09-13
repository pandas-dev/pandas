from __future__ import print_function
from pandas.compat import iteritems

"""
Delayed (rename?)

delayed selection api through magic `X` variable
"""

# includes large portions of pandas_ply, see LICENSES

_error_doc = """
pandas `X` is a deferred object that cannot be passed into
functions. {case} which is invalid.
To pass a deferred Series into a function, use the .pipe
function, for example, X.a.pipe(np.log), instead np.log(df.a)
"""

_disallow_attr = [
    # numpy / pandas access
    '__array_struct__', '__array_interface__',
    '_typ', '_data', 'index', 'columns', 'values',
    # python methods
    '__iter__', '__len__'
]


class Expression(object):
    """
    Expression is the (abstract) base class for symbolic expressions.
    Symbolic expressions are encoded representations of Python expressions,
    kept on ice until you are ready to evaluate them.

    If an expression is complete, it will act as a 1-argument
    function, taking a DataFrame whose context to evaluate
    the expr in.  If not complete, __call__ will create
    a symbolic call node.
    """

    def _eval(self, context, **options):
        """Evaluate a symbolic expression.

        Args:
            context: The context object for evaluation. Currently, this is a
                dictionary mapping symbol names to values,
            `**options`: Options for evaluation. Currently, the only option is
                `log`, which results in some debug output during evaluation if
                it is set to `True`.

        Returns:
            anything
        """
        raise NotImplementedError

    def __repr__(self):
        raise NotImplementedError

    def __getattr__(self, name):
        """Construct a symbolic representation of `getattr(self, name)`."""
        if name in _disallow_attr:
            msg = "The {0} attribuate was called on the object".format(name)
            raise TypeError(_error_doc.format(case=msg))

        # XXX currently assuming complete-ness
        # always alternates, which often works but
        # not always, eg., `X.c.str.upper()` fails
        complete = not self._complete
        return GetAttr(self, name, complete)

    def __getitem__(self, name):
        return GetItem(self, name, complete=True)

    def __call__(self, *args, **kwargs):
        if self._complete:
            # selection lambda passed to pandas
            df, = args
            return self._eval({0: df})
        # symbolic call
        return Call(self, args=args, kwargs=kwargs)

    # error trapping
    def __array__(self, *args, **kwargs):
        msg = "The object was attempted to be converted to a numpy array"
        raise TypeError(_error_doc.format(case=msg))


_allowed_magic_methods = [
    '__abs__', '__and__', '__add__', '__div__', '__radd__', '__div__',
    '__eq__', '__ge__', '__gt__', '__lt__', '__mul__', '__ne__', '__neg__',
    '__or__', '__pos__', '__pow__', '__sub__', '__truediv__'
]


def _get_sym_magic_method(name):
    def magic_method(self, *args, **kwargs):
        return Call(GetAttr(self, name), args, kwargs)
    return magic_method


for name in _allowed_magic_methods:
    setattr(Expression, name, _get_sym_magic_method(name))


class Symbol(Expression):
    """`Symbol(name)` is an atomic symbolic expression, labelled with an
    arbitrary `name`."""

    def __init__(self, name, complete=False):
        self._name = name
        self._complete = complete

    def _eval(self, context, **options):
        if options.get('log'):
            print('Symbol._eval', repr(self))
        result = context[self._name]
        if options.get('log'):
            print('Returning', repr(self), '=>', repr(result))
        return result

    def __repr__(self):
        return 'Symbol(%s)' % repr(self._name)


class GetAttr(Expression):
    """`GetItem(obj, name)` is a symbolic expression representing the result of
    `getattr(obj, name)`. (`obj` and `name` can themselves be symbolic.)"""

    def __init__(self, obj, name, complete=True):
        self._obj = obj
        self._name = name
        self._complete = complete

    def _eval(self, context, **options):
        if options.get('log'):
            print('GetAttr._eval', repr(self))
        evaled_obj = eval_if_symbolic(self._obj, context, **options)
        result = getattr(evaled_obj, self._name)
        if options.get('log'):
            print('Returning', repr(self), '=>', repr(result))
        return result

    def __repr__(self):
        return 'getattr(%s, %s)' % (repr(self._obj), repr(self._name))


class GetItem(Expression):
    """`GetAttr(obj, name)` is a symbolic expression representing the result of
    `getattr(obj, name)`. (`obj` and `name` can themselves be symbolic.)"""

    def __init__(self, obj, name, complete=True):
        self._obj = obj
        self._name = name
        self._complete = complete

    def _eval(self, context, **options):
        if options.get('log'):
            print('GetItem._eval', repr(self))
        evaled_obj = eval_if_symbolic(self._obj, context, **options)
        result = evaled_obj[self._name]
        if options.get('log'):
            print('Returning', repr(self), '=>', repr(result))
        return result

    def __repr__(self):
        return 'getitem(%s, %s)' % (repr(self._obj), repr(self._name))


class Call(Expression):
    """`Call(func, args, kwargs)` is a symbolic expression representing the
    result of `func(*args, **kwargs)`. (`func`, each member of the `args`
    iterable, and each value in the `kwargs` dictionary can themselves be
    symbolic)."""

    def __init__(self, func, args=None, kwargs=None, complete=True):
        self._func = func
        if not args:
            args = []
        if not kwargs:
            kwargs = {}
        self._args = args
        self._kwargs = kwargs
        self._complete = True

    def _eval(self, context, **options):
        if options.get('log'):
            print('Call._eval', repr(self))
        evaled_func = eval_if_symbolic(self._func, context, **options)
        evaled_args = [eval_if_symbolic(v, context, **options)
                       for v in self._args]
        evaled_kwargs = dict((k, eval_if_symbolic(v, context, **options))
                             for k, v in iteritems(self._kwargs))
        result = evaled_func(*evaled_args, **evaled_kwargs)
        if options.get('log'):
            print('Returning', repr(self), '=>', repr(result))
        return result

    def __repr__(self):
        return '{func}(*{args}, **{kwargs})'.format(
            func=repr(self._func),
            args=repr(self._args),
            kwargs=repr(self._kwargs))


def eval_if_symbolic(obj, context, **options):
    """Evaluate an object if it is a symbolic expression, or otherwise just
    returns it back.

    Args:
        obj: Either a symbolic expression, or anything else (in which case this
            is a noop).
        context: Passed as an argument to `obj._eval` if `obj` is symbolic.
        `**options`: Passed as arguments to `obj._eval` if `obj` is symbolic.

    Returns:
        anything

    Examples:
        >>> eval_if_symbolic(Symbol('x'), {'x': 10})
        10
        >>> eval_if_symbolic(7, {'x': 10})
        7
    """
    return obj._eval(context, **options) if hasattr(obj, '_eval') else obj


def to_callable(obj):
    """Turn an object into a callable.

    Args:
        obj: This can be

            * **a symbolic expression**, in which case the output callable
              evaluates the expression with symbols taking values from the
              callable's arguments (listed arguments named according to their
              numerical index, keyword arguments named according to their
              string keys),
            * **a callable**, in which case the output callable is just the
              input object, or
            * **anything else**, in which case the output callable is a
              constant function which always returns the input object.

    Returns:
        callable

    Examples:
        >>> to_callable(Symbol(0) + Symbol('x'))(3, x=4)
        7
        >>> to_callable(lambda x: x + 1)(10)
        11
        >>> to_callable(12)(3, x=4)
        12
    """
    if hasattr(obj, '_eval'):
        return lambda *args, **kwargs: obj._eval(
            dict(enumerate(args), **kwargs))
    elif callable(obj):
        return obj
    else:
        return lambda *args, **kwargs: obj


# keep?
def sym_call(func, *args, **kwargs):
    """Construct a symbolic representation of `func(*args, **kwargs)`.

    This is necessary because `func(symbolic)` will not (ordinarily) know to
    construct a symbolic expression when it receives the symbolic
    expression `symbolic` as a parameter (if `func` is not itself symbolic).
    So instead, we write `sym_call(func, symbolic)`.

    Tip: If the main argument of the function is a (symbolic) DataFrame, then
    pandas' `pipe` method takes care of this problem without `sym_call`. For
    instance, while `np.sqrt(X)` won't work, `X.pipe(np.sqrt)` will.

    Args:
      func: Function to call on evaluation (can be symbolic).
      `*args`: Arguments to provide to `func` on evaluation (can be symbolic).
      `**kwargs`: Keyword arguments to provide to `func` on evaluation (can be
          symbolic).

    Returns:
        `ply.symbolic.Expression`

    Example:
        >>> sym_call(math.sqrt, Symbol('x'))._eval({'x': 16})
        4
    """

    return Call(func, args=args, kwargs=kwargs)

X = Symbol(0)
"""A Symbol for "the first argument" (for convenience)."""
