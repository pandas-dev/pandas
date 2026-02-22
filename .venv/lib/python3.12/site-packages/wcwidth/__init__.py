"""
Wcwidth module.

https://github.com/jquast/wcwidth
"""
# re-export all functions & definitions, even private ones, from top-level
# module path, to allow for 'from wcwidth import _private_func'.  Of course,
# user beware that any _private functions or variables not exported by __all__
# may disappear or change signature at any future version.

# local
from .wcwidth import ZERO_WIDTH  # noqa
from .wcwidth import (WIDE_EASTASIAN,
                      AMBIGUOUS_EASTASIAN,
                      VS16_NARROW_TO_WIDE,
                      clip,
                      ljust,
                      rjust,
                      width,
                      center,
                      wcwidth,
                      wcswidth,
                      list_versions,
                      iter_sequences,
                      strip_sequences,
                      _wcmatch_version,
                      _wcversion_value)
from .bisearch import bisearch as _bisearch
from .grapheme import grapheme_boundary_before  # noqa
from .grapheme import iter_graphemes, iter_graphemes_reverse
from .textwrap import SequenceTextWrapper, wrap
from .sgr_state import propagate_sgr

# The __all__ attribute defines the items exported from statement,
# 'from wcwidth import *', but also to say, "This is the public API".
__all__ = ('wcwidth', 'wcswidth', 'width', 'iter_sequences', 'iter_graphemes',
           'iter_graphemes_reverse', 'grapheme_boundary_before',
           'ljust', 'rjust', 'center', 'wrap', 'clip', 'strip_sequences',
           'list_versions', 'propagate_sgr')

# Using 'hatchling', it does not seem to provide the pyproject.toml nicety, "dynamic = ['version']"
# like flit_core, maybe there is some better way but for now we have to duplicate it in both places
__version__ = '0.6.0'
