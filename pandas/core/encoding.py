#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Mutable sequences handling? specificaly tuples.
# generators - via wrapper?

import unittest

from pandas.util import py3compat
from pandas.core.common import  PandasError
import numpy as np
import sys


try:
    next
except NameError:  # pragma: no cover
    # Python < 2.6
    def next(x):
        return x.next()

# this should live in some package-wide conf object
input_encoding='utf-8'
perform_conversion=True
guess_enc_on_decode_failure=True
guess_enc_min_char_count=100
guess_enc_max_iter=5000
guess_enc_min_confidence=0.8

def set_input_encoding(encoding):
    global input_encoding
    input_encoding=encoding

def _should_process(obj):
    """
    A predicate function which determines whether obj should
    be processed for byte-string conversion based on it's type.

    Parameters
    ----------
    obj - any object

    Returns
    -------
    bool - True if the object should be processed
    """

    # pd.Index* are isinstance(np.ndarray) but should be excluded
    # because their constructors call decode directly.
    #
    return perform_conversion and \
      ( isinstance(obj,(list,dict)) or \
        type(obj) == np.ndarray or \
        type(obj) == np.void )

def _can_import(name):
    """
    Returns True if the named module/package can be imported"

    Parameters
    ----------
    `name` -  package / module name.

    Returns
    -------
    bool - True if `name` can be imported.

    """
    try:
        __import__(name)
        return True
    except ImportError:
        return False

def _decode_obj(obj, encoding):
    """
    Recieves an object, and decodes any non-ascii byte-strings found
    to unicode using the given encoding.

    You should use `decode_catch_errors` to get friendly error messages
    when decoding fails.

    supports str,unicode, and mutable sequences
    This function iterates over `seq`, decoding any byte-string object found using the
    given `encoding`.

    supports str/unicode and mutable sequences as input, all others
    are returned as-is (including generators, for now)

    Handles arbitrarily nested sequences.

    Parameters
    ----------
    `obj` - any object.

    Returns
    -------
    result -  seq with all non-ascii bytestring decoded into utf-8

    Raises
    ------
    UnicodeDecodeError - if decoding with the given encoding fails
    """

    import types
    def _dec_str(s,encoding=encoding):
        try:
            s.encode('ascii') # if it's ascii leave it alone
        except UnicodeDecodeError:
            s = s.decode(encoding) # if not, convert to unicode
            # might raise another UnicodeDecodeError - handled by the caller
        return s

    def _dec_seq(seq):
        if isinstance(seq, dict):
           for k in seq.keys():  # grab the list of keys before we do mutation
               v=seq[k]
               if isinstance(k, str):
                   k = _dec_str(k)
               elif _should_process(k): # keys are immutable, need this?
                   k = (yield _dec_seq(k))

               if isinstance(v, str):
                   v = _dec_str(v)
               elif _should_process(v):
                   v = (yield _dec_seq(v))

               seq.pop(k) # need this
               seq[k] = v

        else:
            for i,e in enumerate(seq):
                if isinstance(e, str):
                    seq[i] = _dec_str(e)
                elif _should_process(e):
                    (yield _dec_seq(e))

        yield seq

    if py3compat.PY3:
        return obj

    if isinstance(obj, basestring): # strings are simple
        if isinstance(obj, str):
            obj=_dec_str(obj)
        return obj

    if not _should_process(obj): # misc. objects are too
        return obj

    s = [_dec_seq(obj)]
    values = []
    while True: # others - not so much, let's see what we can do.
        g = s.pop()
        if values:
            e = g.send(values.pop())
        else:
            e = next(g)
        if type(e) == types.GeneratorType:
            s.extend([g, e])
        else:
            if s:
                values.append(e)
            else:
                return e

def _extract_txt_from_obj(obj,max_iter=sys.maxint):
    """
    a generator which walks `obj`, yielding any byte-string found

    will stop after at most `max_iter` iterations.

    Parameters
    ----------
    `obj` - any iterable

    Yields
    -------
    byte-strings.

    Raises
    ------
    StopIteration - when there are no more byte-strings in the sequence

    """

    if obj is None or isinstance(obj,basestring):
        if isinstance(obj,unicode):
            return
        yield obj
        return

    s = [iter(obj)]
    cnt=0
    while s:
        g = s.pop()
        for e in g:
            cnt+=1
            if isinstance(e, str):
                yield e
            elif isinstance(e, dict):
                s.extend([g, e.iterkeys(), e.itervalues()])
            elif _should_process(e):
                s.append(g, iter(e))

            if cnt >= max_iter:
                return

def _detect_encoding(obj,min_cnt=guess_enc_min_char_count,max_iter=guess_enc_max_iter):
    """
    extracts byte-string from obj via `_extract_txt_from_obj` and uses
    the `chardet` package to detect the encoding used.

    Can handle nested sequences, also looks at dict keys and values.

    Parameters
    ----------
    `obj` - input string or sequence

    `min_cnt` - specifies the minimum amount of characters which must be fed to to
    the detector before we allow a decision.

    `max_iter` - an upper bound on the number of elements examined in the sequence
    looking for text.
    This guards against the corner-case of a huge list with a decoding error only
    near it's end.

    Returns
    -------
    `result` - {'encoding': str, 'confidence': float}, or
             {} if no encoding was found.
    """
    if not _can_import("chardet"):
        return {}

    from chardet.universaldetector import UniversalDetector
    detector = UniversalDetector()
    cnt = 0 # keep track of number of characters processed
    for txt in _extract_txt_from_obj(obj,max_iter):
        cnt += len(txt)
        detector.feed(txt)
        if (cnt > min_cnt and detector.done) :
            break
    detector.close()
    res=detector.result
    if res and res['confidence'] > guess_enc_min_confidence\
      and cnt > min_cnt:
        return detector.result
    else:
        return {}

def decode_catch_errors(obj, encoding=None):
    """
    Delegates to `_decode_obj` in order to convert byte-string within obj into
    unicode when necessary. If a decode error occurs, prints a user friendly
    error message, and if the chardet library is available will try to give
    the user a good guesss about the encoding used by extracting text from `obj`

    Parameters
    ----------
    `obj` - anything
    encoding - an acceptable encoding to be passed to str.decode()

    Returns
    -------
    `result` - `obj` with byte-strings decoded into unicode strings

    Raises
    ------
    `PandasError(msg)` - with msg being a friendly error message to the user
    """

    try:
        encoding = encoding or input_encoding
        return _decode_obj(obj, encoding)
    except UnicodeDecodeError:
        from textwrap import dedent
        msg = \
            """
        The input Data contains strings that cannot be decoded with `%s`.
        You should specify a correct encoding to the object constructor,
        or set the value of the default input encoding in XXX.
        """

        s = dedent(msg) % encoding
        if guess_enc_on_decode_failure:
            if not _can_import("chardet"):
                s += 'The "chardet" package is not installed - ' +\
                     "can't suggest an encoding."
            else:
                det_enc=_detect_encoding(obj)
                if det_enc:
                    conf = det_enc['confidence']
                    enc = det_enc['encoding']
                    s += 'You might try "%s" as the encoding (Confidence: %2.1f)'\
                  % (enc, conf)

        raise PandasError(s)
