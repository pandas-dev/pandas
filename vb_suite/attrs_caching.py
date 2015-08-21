from vbench.benchmark import Benchmark

common_setup = """from .pandas_vb_common import *
"""

#----------------------------------------------------------------------
# DataFrame.index / columns property lookup time

setup = common_setup + """
df = DataFrame(np.random.randn(10, 6))
cur_index = df.index
"""
stmt = "foo = df.index"

getattr_dataframe_index = Benchmark(stmt, setup,
                                    name="getattr_dataframe_index")

stmt = "df.index = cur_index"
setattr_dataframe_index = Benchmark(stmt, setup,
                                    name="setattr_dataframe_index")
