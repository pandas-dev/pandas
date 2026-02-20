__all__ = []

from .ndindex import ndindex

__all__ += ['ndindex']

from .shapetools import broadcast_shapes, iter_indices, AxisError, BroadcastError

__all__ += ['broadcast_shapes', 'iter_indices', 'AxisError', 'BroadcastError']

from .slice import Slice

__all__ += ['Slice']

from .integer import Integer

__all__ += ['Integer']

from .tuple import Tuple

__all__ += ['Tuple']

from .ellipsis import ellipsis

__all__ += ['ellipsis']

from .newaxis import Newaxis

__all__ += ['Newaxis']

from .integerarray import IntegerArray

__all__ += ['IntegerArray']

from .booleanarray import BooleanArray

__all__ += ['BooleanArray']

from .chunking import ChunkSize

__all__ += ['ChunkSize']

from ._version import get_versions
__version__ = get_versions()['version']
del get_versions
