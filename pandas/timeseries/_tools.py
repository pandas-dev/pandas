"""
A collection of helper tools


Decorators
==========

.. autoclass:: docwrapper

.. autoclass:: deprecated_for

"""

import functools
import warnings



class docwrapper(object):
    """
    Decorator that updates the docstring of a function with a dictionary
    of templates.

    Parameters
    ----------
    template : dictionary
        A dictionary giving for each key the docstring to extend.
    """
    #
    def __init__(self, template):
        self.template = template
    #
    def __call__(self, func):
        def wrapped(*args, **kwargs):
            "Just call the function w/ the proper arguments"
            return func(*args, **kwargs)
        wrapped.__name__ = func.__name__
        wrapped.__dict__ = func.__dict__
        wrapped.__doc__ = ((func.__doc__ or "") % self.template) or None
        return wrapped




class deprecated_for:
    """
    Decorator marking a function as deprecated.

    When a decorated function is called, a warning is emitted.

    Parameters
    ----------
    newfunc : function, optional
        Function that should replace the deprecated one.

    Examples
    --------
    >>> def new_function(*args, **kwargs):
    ...     do_something
    ...
    >>> @deprecated_for(new_function)
    >>> def old_function(*args, **kwargs):
    ...    do_something
    ...
    >>> old_function(*args)
    DeprecationWarning: The function `old_function` is deprecated.
    Please use the `new_function` instead.
    """
    #
    def __init__(self, newfunc=None):
        self.replacement = newfunc
        msg = "The function `%(oldname)s` is deprecated.\n"
        doctemplate = "(Deprecated function)\n%(olddoc)s"
        if newfunc is not None:
            msg += "Please use the `%s` function instead." % newfunc.__name__
            newdoc = newfunc.__doc__
            if newdoc is not None:
                doctemplate += "\n(New usage)\n%s" % newdoc
        self.msg = msg
        self.doctemplate = doctemplate
    #
    def __call__(self, func):
        oldinfo = dict(oldname=func.__name__, olddoc=func.__doc__ or "")
        msg = self.msg % oldinfo
        def wrapped(*args, **kwargs):
            warnings.warn_explicit(msg, category=DeprecationWarning,
                                   filename=func.func_code.co_filename,
                                   lineno=func.func_code.co_firstlineno + 1)
            return func(*args, **kwargs)
        wrapped.__name__ = func.__name__
        wrapped.__dict__ = func.__dict__
        wrapped.__doc__ = (self.doctemplate % oldinfo) or None
        return wrapped

deprecated = deprecated_for()
