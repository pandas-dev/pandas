"""Module for scope operations
"""

import sys
import operator
import struct
import inspect
import datetime
import itertools
import pprint

import pandas as pd
from pandas.compat import DeepChainMap, map
from pandas.core import common as com
from pandas.core.base import StringMixin
from pandas.computation.ops import UndefinedVariableError


def _ensure_scope(level, global_dict=None, local_dict=None, resolvers=(),
                  target=None, **kwargs):
    """Ensure that we are grabbing the correct scope."""
    return Scope(level + 1, gbls=global_dict, lcls=local_dict,
                 resolvers=resolvers, target=target)


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



_DEFAULT_GLOBALS = {
    'Timestamp': pd.lib.Timestamp,
    'datetime': datetime.datetime,
    'True': True,
    'False': False,
    'list': list,
    'tuple': tuple
}


def _is_resolver(x):
    return isinstance(x, Resolver)


def _get_pretty_string(obj):
    sio = StringIO()
    pprint.pprint(obj, stream=sio)
    return sio.getvalue()


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
    __slots__ = 'level', 'scope', 'target', 'ntemps'

    def __init__(self, level, gbls=None, lcls=None, resolvers=(), target=None):
        self.level = level + 1

        # shallow copy because we don't want to keep filling this up with what
        # was there before if there are multiple calls to Scope/_ensure_scope
        self.scope = DeepChainMap(_DEFAULT_GLOBALS.copy())
        self.target = target
        self.ntemps = 0  # number of temporary variables in this scope

        if isinstance(lcls, Scope):
            self.scope.update(lcls.scope)
            if lcls.target is not None:
                self.target = lcls.target
            self.update(lcls.level)

        frame = sys._getframe(self.level)

        try:
            # shallow copy here because we don't want to replace what's in
            # scope when we align terms (alignment accesses the underlying
            # numpy array of pandas objects)
            if not isinstance(lcls, Scope):
                self.scope = self.scope.new_child((lcls or frame.f_locals).copy())
            self.scope = self.scope.new_child((gbls or frame.f_globals).copy())
        finally:
            del frame

        # assumes that resolvers are going from outermost scope to inner
        if isinstance(lcls, Scope):
            resolvers += tuple(lcls.resolvers.maps)
        self.resolvers = DeepChainMap(*resolvers)

    def __unicode__(self):
        scope_keys = _get_pretty_string(self.scope.keys())
        res_keys = _get_pretty_string(self.resolvers.keys())
        return 'Scope(scope=%s, resolvers=%s)' % (scope_keys, res_keys)

    @property
    def has_resolvers(self):
        return bool(self.nresolvers)

    @property
    def nresolvers(self):
        return len(self.resolvers)

    def resolve(self, key, is_local):
        """Resolve a variable name in a possibly local context

        Parameters
        ----------
        key : text_type
            A variable name
        is_local : bool
            Flag indicating whether the variable is local or not (prefixed with
            the '@' symbol)

        Returns
        -------
        value : object
            The value of a particular variable
        """
        try:
            # only look for locals in outer scope
            if is_local:
                return self.scope[key]

            # not a local variable so check in resolvers if we have them
            if self.has_resolvers:
                return self.resolvers[key]

            # if we're here that means that we have no locals and we also have
            # no resolvers
            assert not is_local and not self.has_resolvers
            return self.scope[key]
        except KeyError:
            raise UndefinedVariableError(key)

    def swapkey(self, old_key, new_key, new_value=None):
        if self.has_resolvers:
            maps = self.resolvers.maps + self.scope.maps
        else:
            maps = self.scope.maps

        for mapping in maps:
            if old_key in mapping:
                if new_value is None:
                    mapping[new_key] = mapping.pop(old_key)
                else:
                    mapping[new_key] = new_value
                return
        raise KeyError(old_key)

    def _get_vars(self, stack, scopes):
        variables = itertools.product(scopes, stack)
        for scope, (frame, _, _, _, _, _) in variables:
            try:
                d = getattr(frame, 'f_' + scope)
                self.scope = self.scope.new_child(d)
            finally:
                # won't remove it, but DECREF it
                # in Py3 this probably isn't necessary since frame won't be
                # scope after the loop
                del frame

    def update(self, level):
        """Update the current scope by going back `level` levels.

        Parameters
        ----------
        level : int or None, optional, default None
        """
        sl = level + 1

        # add sl frames to the scope starting with the
        # most distant and overwriting with more current
        # makes sure that we can capture variable scope
        stack = inspect.stack()

        try:
            self._get_vars(stack[:sl], scopes=['locals'])
        finally:
            del stack[:], stack

    def add_tmp(self, value):
        """Add a temporary variable to the scope.

        Parameters
        ----------
        value : object
            An arbitrary object to be assigned to a temporary variable.

        Returns
        -------
        name : basestring
            The name of the temporary variable created.
        """
        name = 'tmp_var_{0}_{1}_{2}'.format(type(value).__name__, self.ntemps,
                                            _raw_hex_id(self))

        # add to inner most scope
        assert name not in self.scope.maps[0]
        self.scope.maps[0][name] = value

        # only increment if the variable gets put in the scope
        self.ntemps += 1
        return name

    def remove_tmp(self, name):
        del self.scope[name]
        self.ntemps -= 1
