import inspect
import logging
import sys
import threading
from functools import wraps
from itertools import chain

from numba.core import config


class TLS(threading.local):
    """Use a subclass to properly initialize the TLS variables in all threads.""" # noqa: E501

    def __init__(self):
        self.tracing = False
        self.indent = 0


tls = TLS()


def find_function_info(func, spec, args):
    """Return function meta-data in a tuple.

    (name, type)"""

    module = getattr(func, "__module__", None)
    name = getattr(func, "__name__", None)
    self = getattr(func, "__self__", None)
    cname = None
    if self:
        cname = self.__name__
        # cname = self.__class__.__name__
    # Try to deduce the class' name even for unbound methods from their
    # first argument, which we assume to be a class instance if named 'self'...
    elif len(spec.args) and spec.args[0] == "self":
        cname = args[0].__class__.__name__
    # ...or a class object if named 'cls'
    elif len(spec.args) and spec.args[0] == "cls":
        cname = args[0].__name__
    if name:
        qname = []
        if module and module != "__main__":
            qname.append(module)
            qname.append(".")
        if cname:
            qname.append(cname)
            qname.append(".")
        qname.append(name)
        name = "".join(qname)
    return name, None


def chop(value):
    MAX_SIZE = 320
    s = repr(value)
    if len(s) > MAX_SIZE:
        return s[:MAX_SIZE] + "..." + s[-1]
    else:
        return s


def create_events(fname, spec, args, kwds):
    values = dict()
    if spec.defaults:
        values = dict(zip(spec.args[-len(spec.defaults) :], spec.defaults))
    values.update(kwds)
    values.update(list(zip(spec.args[: len(args)], args)))
    positional = ["%s=%r" % (a, values.pop(a)) for a in spec.args]
    anonymous = [str(a) for a in args[len(positional) :]]
    keywords = ["%s=%r" % (k, values[k]) for k in sorted(values.keys())]
    params = ", ".join([f for f in chain(positional, anonymous, keywords) if f])

    enter = [">> ", tls.indent * " ", fname, "(", params, ")"]
    leave = ["<< ", tls.indent * " ", fname]
    return enter, leave


def dotrace(*args, **kwds):
    """Function decorator to trace a function's entry and exit.

    *args: categories in which to trace this function. Example usage:

    @trace
    def function(...):...

    @trace('mycategory')
    def function(...):...


    """

    recursive = kwds.get("recursive", False)

    def decorator(func):
        spec = None
        logger = logging.getLogger("trace")

        def wrapper(*args, **kwds):
            if not logger.isEnabledFor(logging.INFO) or tls.tracing:
                return func(*args, **kwds)

            fname, ftype = find_function_info(func, spec, args)

            try:
                tls.tracing = True
                enter, leave = create_events(fname, spec, args, kwds)

                try:
                    logger.info("".join(enter))
                    tls.indent += 1
                    try:
                        try:
                            tls.tracing = False
                            result = func(*args, **kwds)
                        finally:
                            tls.tracing = True
                    except: # noqa: E722
                        type, value, traceback = sys.exc_info()
                        leave.append(" => exception thrown\n\traise ")
                        mname = type.__module__
                        if mname != "__main__":
                            leave.append(mname)
                            leave.append(".")
                        leave.append(type.__name__)
                        if value.args:
                            leave.append("(")
                            leave.append(", ".join(chop(v) for v in value.args))
                            leave.append(")")
                        else:
                            leave.append("()")
                        raise
                    else:
                        if result is not None:
                            leave.append(" -> ")
                            leave.append(chop(result))
                finally:
                    tls.indent -= 1
                    logger.info("".join(leave))
            finally:
                tls.tracing = False
            return result

        # wrapper end

        rewrap = lambda x: x
        # Unwrap already wrapped functions
        # (to be rewrapped again later)
        if isinstance(func, classmethod):
            rewrap = type(func)
            # Note: 'func.__func__' only works in Python 3
            func = func.__get__(True).__func__
        elif isinstance(func, staticmethod):
            rewrap = type(func)
            # Note: 'func.__func__' only works in Python 3
            func = func.__get__(True)
        elif isinstance(func, property):
            raise NotImplementedError

        spec = inspect.getfullargspec(func)
        return rewrap(wraps(func)(wrapper))

    arg0 = len(args) and args[0] or None
    # not supported yet...
    if recursive:
        raise NotImplementedError
        if inspect.ismodule(arg0):
            for n, f in inspect.getmembers(arg0, inspect.isfunction):
                setattr(arg0, n, decorator(f))
            for n, c in inspect.getmembers(arg0, inspect.isclass):
                dotrace(c, *args, recursive=recursive)
        elif inspect.isclass(arg0):
            for n, f in inspect.getmembers(
                arg0, lambda x: (inspect.isfunction(x) or inspect.ismethod(x))
            ):
                setattr(arg0, n, decorator(f))

    if callable(arg0) or type(arg0) in (classmethod, staticmethod):
        return decorator(arg0)
    elif isinstance(arg0, property):
        # properties combine up to three functions: 'get', 'set', 'del',
        # so let's wrap them all.
        pget, pset, pdel = None, None, None
        if arg0.fget:
            pget = decorator(arg0.fget)
        if arg0.fset:
            pset = decorator(arg0.fset)
        if arg0.fdel:
            pdel = decorator(arg0.fdel)
        return property(pget, pset, pdel)

    else:
        return decorator


def notrace(*args, **kwds):
    """Just a no-op in case tracing is disabled."""

    def decorator(func):
        return func

    arg0 = len(args) and args[0] or None

    if callable(arg0) or type(arg0) in (classmethod, staticmethod):
        return decorator(arg0)
    else:
        return decorator


def doevent(msg):
    msg = ["== ", tls.indent * " ", msg]
    logger = logging.getLogger("trace")
    logger.info("".join(msg))


def noevent(msg):
    pass


if config.TRACE:
    logger = logging.getLogger("trace")
    logger.setLevel(logging.INFO)
    logger.handlers = [logging.StreamHandler()]
    trace = dotrace
    event = doevent
else:
    trace = notrace
    event = noevent
