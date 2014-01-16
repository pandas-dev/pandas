from pandas.compat import StringIO, callable
from pandas.lib import cache_readonly
import sys
import warnings


def deprecate(name, alternative, alt_name=None):
    alt_name = alt_name or alternative.__name__

    def wrapper(*args, **kwargs):
        warnings.warn("%s is deprecated. Use %s instead" % (name, alt_name),
                      FutureWarning)
        return alternative(*args, **kwargs)
    return wrapper

# Substitution and Appender are derived from matplotlib.docstring (1.1.0)
# module http://matplotlib.sourceforge.net/users/license.html


class Substitution(object):
    """
    A decorator to take a function's docstring and perform string
    substitution on it.

    This decorator should be robust even if func.__doc__ is None
    (for example, if -OO was passed to the interpreter)

    Usage: construct a docstring.Substitution with a sequence or
    dictionary suitable for performing substitution; then
    decorate a suitable function with the constructed object. e.g.

    sub_author_name = Substitution(author='Jason')

    @sub_author_name
    def some_function(x):
        "%(author)s wrote this function"

    # note that some_function.__doc__ is now "Jason wrote this function"

    One can also use positional arguments.

    sub_first_last_names = Substitution('Edgar Allen', 'Poe')

    @sub_first_last_names
    def some_function(x):
        "%s %s wrote the Raven"
    """
    def __init__(self, *args, **kwargs):
        if (args and kwargs):
            raise AssertionError( "Only positional or keyword args are allowed")

        self.params = args or kwargs

    def __call__(self, func):
        func.__doc__ = func.__doc__ and func.__doc__ % self.params
        return func

    def update(self, *args, **kwargs):
        "Assume self.params is a dict and update it with supplied args"
        self.params.update(*args, **kwargs)

    @classmethod
    def from_params(cls, params):
        """
        In the case where the params is a mutable sequence (list or dictionary)
        and it may change before this class is called, one may explicitly use a
        reference to the params rather than using *args or **kwargs which will
        copy the values and not reference them.
        """
        result = cls()
        result.params = params
        return result


class Appender(object):
    """
    A function decorator that will append an addendum to the docstring
    of the target function.

    This decorator should be robust even if func.__doc__ is None
    (for example, if -OO was passed to the interpreter).

    Usage: construct a docstring.Appender with a string to be joined to
    the original docstring. An optional 'join' parameter may be supplied
    which will be used to join the docstring and addendum. e.g.

    add_copyright = Appender("Copyright (c) 2009", join='\n')

    @add_copyright
    def my_dog(has='fleas'):
        "This docstring will have a copyright below"
        pass
    """
    def __init__(self, addendum, join='', indents=0):
        if indents > 0:
            self.addendum = indent(addendum, indents=indents)
        else:
            self.addendum = addendum
        self.join = join

    def __call__(self, func):
        func.__doc__ = func.__doc__ if func.__doc__ else ''
        self.addendum = self.addendum if self.addendum else ''
        docitems = [func.__doc__, self.addendum]
        func.__doc__ = self.join.join(docitems)
        return func


def indent(text, indents=1):
    if not text or not isinstance(text, str):
        return ''
    jointext = ''.join(['\n'] + ['    '] * indents)
    return jointext.join(text.split('\n'))


def suppress_stdout(f):
    def wrapped(*args, **kwargs):
        try:
            sys.stdout = StringIO()
            f(*args, **kwargs)
        finally:
            sys.stdout = sys.__stdout__

    return wrapped


class KnownFailureTest(Exception):
    '''Raise this exception to mark a test as a known failing test.'''
    pass


def knownfailureif(fail_condition, msg=None):
    """
    Make function raise KnownFailureTest exception if given condition is true.

    If the condition is a callable, it is used at runtime to dynamically
    make the decision. This is useful for tests that may require costly
    imports, to delay the cost until the test suite is actually executed.

    Parameters
    ----------
    fail_condition : bool or callable
        Flag to determine whether to mark the decorated test as a known
        failure (if True) or not (if False).
    msg : str, optional
        Message to give on raising a KnownFailureTest exception.
        Default is None.

    Returns
    -------
    decorator : function
        Decorator, which, when applied to a function, causes SkipTest
        to be raised when `skip_condition` is True, and the function
        to be called normally otherwise.

    Notes
    -----
    The decorator itself is decorated with the ``nose.tools.make_decorator``
    function in order to transmit function name, and various other metadata.

    """
    if msg is None:
        msg = 'Test skipped due to known failure'

    # Allow for both boolean or callable known failure conditions.
    if callable(fail_condition):
        fail_val = fail_condition
    else:
        fail_val = lambda: fail_condition

    def knownfail_decorator(f):
        # Local import to avoid a hard nose dependency and only incur the
        # import time overhead at actual test-time.
        import nose

        def knownfailer(*args, **kwargs):
            if fail_val():
                raise KnownFailureTest(msg)
            else:
                return f(*args, **kwargs)
        return nose.tools.make_decorator(f)(knownfailer)

    return knownfail_decorator
