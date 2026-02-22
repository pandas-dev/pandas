# -*- coding: utf-8 -*-
#
# python-json-pointer - An implementation of the JSON Pointer syntax
# https://github.com/stefankoegl/python-json-pointer
#
# Copyright (c) 2011 Stefan Kögl <stefan@skoegl.net>
# All rights reserved.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions
# are met:
#
# 1. Redistributions of source code must retain the above copyright
# notice, this list of conditions and the following disclaimer.
# 2. Redistributions in binary form must reproduce the above copyright
# notice, this list of conditions and the following disclaimer in the
# documentation and/or other materials provided with the distribution.
# 3. The name of the author may not be used to endorse or promote products
# derived from this software without specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY THE AUTHOR ``AS IS'' AND ANY EXPRESS OR
# IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES
# OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED.
# IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR ANY DIRECT, INDIRECT,
# INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT
# NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
# DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
# THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
# (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF
# THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
#

""" Identify specific nodes in a JSON document (RFC 6901) """

# Will be parsed by setup.py to determine package metadata
__author__ = 'Stefan Kögl <stefan@skoegl.net>'
__version__ = '3.0.0'
__website__ = 'https://github.com/stefankoegl/python-json-pointer'
__license__ = 'Modified BSD License'

import copy
import re
from collections.abc import Mapping, Sequence
from itertools import tee, chain

_nothing = object()


def set_pointer(doc, pointer, value, inplace=True):
    """Resolves a pointer against doc and sets the value of the target within doc.

    With inplace set to true, doc is modified as long as pointer is not the
    root.

    >>> obj = {'foo': {'anArray': [ {'prop': 44}], 'another prop': {'baz': 'A string' }}}

    >>> set_pointer(obj, '/foo/anArray/0/prop', 55) == \
    {'foo': {'another prop': {'baz': 'A string'}, 'anArray': [{'prop': 55}]}}
    True

    >>> set_pointer(obj, '/foo/yet another prop', 'added prop') == \
    {'foo': {'another prop': {'baz': 'A string'}, 'yet another prop': 'added prop', 'anArray': [{'prop': 55}]}}
    True

    >>> obj = {'foo': {}}
    >>> set_pointer(obj, '/foo/a%20b', 'x') == \
    {'foo': {'a%20b': 'x' }}
    True
    """

    pointer = JsonPointer(pointer)
    return pointer.set(doc, value, inplace)


def resolve_pointer(doc, pointer, default=_nothing):
    """ Resolves pointer against doc and returns the referenced object

    >>> obj = {'foo': {'anArray': [ {'prop': 44}], 'another prop': {'baz': 'A string' }}, 'a%20b': 1, 'c d': 2}

    >>> resolve_pointer(obj, '') == obj
    True

    >>> resolve_pointer(obj, '/foo') == obj['foo']
    True

    >>> resolve_pointer(obj, '/foo/another prop') == obj['foo']['another prop']
    True

    >>> resolve_pointer(obj, '/foo/another prop/baz') == obj['foo']['another prop']['baz']
    True

    >>> resolve_pointer(obj, '/foo/anArray/0') == obj['foo']['anArray'][0]
    True

    >>> resolve_pointer(obj, '/some/path', None) == None
    True

    >>> resolve_pointer(obj, '/a b', None) == None
    True

    >>> resolve_pointer(obj, '/a%20b') == 1
    True

    >>> resolve_pointer(obj, '/c d') == 2
    True

    >>> resolve_pointer(obj, '/c%20d', None) == None
    True
    """

    pointer = JsonPointer(pointer)
    return pointer.resolve(doc, default)


def pairwise(iterable):
    """ Transforms a list to a list of tuples of adjacent items

    s -> (s0,s1), (s1,s2), (s2, s3), ...

    >>> list(pairwise([]))
    []

    >>> list(pairwise([1]))
    []

    >>> list(pairwise([1, 2, 3, 4]))
    [(1, 2), (2, 3), (3, 4)]
    """
    a, b = tee(iterable)
    for _ in b:
        break
    return zip(a, b)


class JsonPointerException(Exception):
    pass


class EndOfList(object):
    """Result of accessing element "-" of a list"""

    def __init__(self, list_):
        self.list_ = list_

    def __repr__(self):
        return '{cls}({lst})'.format(cls=self.__class__.__name__,
                                     lst=repr(self.list_))


class JsonPointer(object):
    """A JSON Pointer that can reference parts of a JSON document"""

    # Array indices must not contain:
    # leading zeros, signs, spaces, decimals, etc
    _RE_ARRAY_INDEX = re.compile('0|[1-9][0-9]*$')
    _RE_INVALID_ESCAPE = re.compile('(~[^01]|~$)')

    def __init__(self, pointer):

        # validate escapes
        invalid_escape = self._RE_INVALID_ESCAPE.search(pointer)
        if invalid_escape:
            raise JsonPointerException('Found invalid escape {}'.format(
                invalid_escape.group()))

        parts = pointer.split('/')
        if parts.pop(0) != '':
            raise JsonPointerException('Location must start with /')

        parts = [unescape(part) for part in parts]
        self.parts = parts

    def to_last(self, doc):
        """Resolves ptr until the last step, returns (sub-doc, last-step)"""

        if not self.parts:
            return doc, None

        for part in self.parts[:-1]:
            doc = self.walk(doc, part)

        return doc, JsonPointer.get_part(doc, self.parts[-1])

    def resolve(self, doc, default=_nothing):
        """Resolves the pointer against doc and returns the referenced object"""

        for part in self.parts:

            try:
                doc = self.walk(doc, part)
            except JsonPointerException:
                if default is _nothing:
                    raise
                else:
                    return default

        return doc

    get = resolve

    def set(self, doc, value, inplace=True):
        """Resolve the pointer against the doc and replace the target with value."""

        if len(self.parts) == 0:
            if inplace:
                raise JsonPointerException('Cannot set root in place')
            return value

        if not inplace:
            doc = copy.deepcopy(doc)

        (parent, part) = self.to_last(doc)

        if isinstance(parent, Sequence) and part == '-':
            parent.append(value)
        else:
            parent[part] = value

        return doc

    @classmethod
    def get_part(cls, doc, part):
        """Returns the next step in the correct type"""

        if isinstance(doc, Mapping):
            return part

        elif isinstance(doc, Sequence):

            if part == '-':
                return part

            if not JsonPointer._RE_ARRAY_INDEX.match(str(part)):
                raise JsonPointerException("'%s' is not a valid sequence index" % part)

            return int(part)

        elif hasattr(doc, '__getitem__'):
            # Allow indexing via ducktyping
            # if the target has defined __getitem__
            return part

        else:
            raise JsonPointerException("Document '%s' does not support indexing, "
                                       "must be mapping/sequence or support __getitem__" % type(doc))

    def get_parts(self):
        """Returns the list of the parts. For example, JsonPointer('/a/b').get_parts() == ['a', 'b']"""

        return self.parts

    def walk(self, doc, part):
        """ Walks one step in doc and returns the referenced part """

        part = JsonPointer.get_part(doc, part)

        assert hasattr(doc, '__getitem__'), "invalid document type %s" % (type(doc),)

        if isinstance(doc, Sequence):
            if part == '-':
                return EndOfList(doc)

            try:
                return doc[part]

            except IndexError:
                raise JsonPointerException("index '%s' is out of bounds" % (part,))

        # Else the object is a mapping or supports __getitem__(so assume custom indexing)
        try:
            return doc[part]

        except KeyError:
            raise JsonPointerException("member '%s' not found in %s" % (part, doc))

    def contains(self, ptr):
        """ Returns True if self contains the given ptr """
        return self.parts[:len(ptr.parts)] == ptr.parts

    def __contains__(self, item):
        """ Returns True if self contains the given ptr """
        return self.contains(item)

    def join(self, suffix):
        """ Returns a new JsonPointer with the given suffix append to this ptr """
        if isinstance(suffix, JsonPointer):
            suffix_parts = suffix.parts
        elif isinstance(suffix, str):
            suffix_parts = JsonPointer(suffix).parts
        else:
            suffix_parts = suffix
        try:
            return JsonPointer.from_parts(chain(self.parts, suffix_parts))
        except:  # noqa E722
            raise JsonPointerException("Invalid suffix")

    def __truediv__(self, suffix):  # Python 3
        return self.join(suffix)

    @property
    def path(self):
        """Returns the string representation of the pointer

        >>> ptr = JsonPointer('/~0/0/~1').path == '/~0/0/~1'
        """
        parts = [escape(part) for part in self.parts]
        return ''.join('/' + part for part in parts)

    def __eq__(self, other):
        """Compares a pointer to another object

        Pointers can be compared by comparing their strings (or splitted
        strings), because no two different parts can point to the same
        structure in an object (eg no different number representations)
        """

        if not isinstance(other, JsonPointer):
            return False

        return self.parts == other.parts

    def __hash__(self):
        return hash(tuple(self.parts))

    def __str__(self):
        return self.path

    def __repr__(self):
        return type(self).__name__ + "(" + repr(self.path) + ")"

    @classmethod
    def from_parts(cls, parts):
        """Constructs a JsonPointer from a list of (unescaped) paths

        >>> JsonPointer.from_parts(['a', '~', '/', 0]).path == '/a/~0/~1/0'
        True
        """
        parts = [escape(str(part)) for part in parts]
        ptr = cls(''.join('/' + part for part in parts))
        return ptr


def escape(s):
    return s.replace('~', '~0').replace('/', '~1')


def unescape(s):
    return s.replace('~1', '/').replace('~0', '~')
