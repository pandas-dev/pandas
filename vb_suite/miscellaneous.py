from vbench.benchmark import Benchmark
from datetime import datetime

common_setup = """from pandas_vb_common import *
"""

#----------------------------------------------------------------------
# cache_readonly

setup = common_setup + """
from pandas.util.decorators import cache_readonly

class Foo:

    @cache_readonly
    def prop(self):
        return 5
obj = Foo()
"""
misc_cache_readonly = Benchmark("obj.prop", setup, name="misc_cache_readonly",
                                ncalls=2000000)

#----------------------------------------------------------------------
# match

setup = common_setup + """
from pandas.util.testing import rands

uniques = np.array([rands(10) for _ in xrange(1000)], dtype='O')
all = uniques.repeat(10)
"""

match_strings = Benchmark("match(all, uniques)", setup,
                          start_date=datetime(2012, 5, 12))
