"""
This module defines the necessary classes and type hints for exposing the bare
minimum of the internal implementations for the pattern (regular expression)
matching backends. The exact structure of the backends is not solidified and is
subject to change.
"""

from collections.abc import (
	Sequence)
from typing import (
	Callable,  # Replaced by `collections.abc.Callable` in 3.9.2.
	Literal,
	Optional)  # Replaced by `X | None` in 3.10.

from .pattern import (
	Pattern)

BackendNamesHint = Literal['best', 'hyperscan', 're2', 'simple']
"""
The supported backend values.
"""

_TestBackendFactoryHint = Optional[Callable[[Sequence[Pattern]], '_Backend']]
"""
Type hint for the test backend factory argument.
"""


class _Backend(object):
	"""
	.. warning:: This class is not part of the public API. It is subject to
		change.

	The :class:`_Backend` class is the abstract base class defining how to match
	files against patterns.
	"""

	def match_file(self, file: str) -> tuple[Optional[bool], Optional[int]]:
		"""
		Check the file against the patterns.

		*file* (:class:`str`) is the normalized file path to check.

		Returns a :class:`tuple` containing whether to include *file* (:class:`bool`
		or :data:`None`), and the index of the last matched pattern (:class:`int` or
		:data:`None`).
		"""
		raise NotImplementedError((
			f"{self.__class__.__module__}.{self.__class__.__qualname__}.match_file() "
			f"must be implemented."
		))  # NotImplementedError
