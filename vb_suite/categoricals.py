from vbench.benchmark import Benchmark
from datetime import datetime

common_setup = """from .pandas_vb_common import *
"""

#----------------------------------------------------------------------
# Series constructors

setup = common_setup + """
s = pd.Series(list('aabbcd') * 1000000).astype('category')
"""

concat_categorical = \
    Benchmark("concat([s, s])", setup=setup, name='concat_categorical',
              start_date=datetime(year=2015, month=7, day=15))
