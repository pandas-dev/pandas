from bottleneck.benchmark.bench import bench
from bottleneck.benchmark.bench_detailed import bench_detailed
from bottleneck.tests.util import get_functions

from . import slow
from ._pytesttester import PytestTester
from .move import (move_argmax, move_argmin, move_max, move_mean, move_median,
                   move_min, move_rank, move_std, move_sum, move_var)
from .nonreduce import replace
from .nonreduce_axis import (argpartition, nanrankdata, partition, push,
                             rankdata)
from .reduce import (allnan, anynan, median, nanargmax, nanargmin, nanmax,
                     nanmean, nanmedian, nanmin, nanstd, nansum, nanvar, ss)

test = PytestTester(__name__)
del PytestTester

from ._version import get_versions  # noqa: E402

__version__ = get_versions()["version"]
del get_versions

from . import _version
__version__ = _version.get_versions()['version']
