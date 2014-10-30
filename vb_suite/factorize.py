from vbench.api import Benchmark
from datetime import datetime

START_DATE = datetime(2014, 10, 13)

# GH 8524

common_setup = """from pandas_vb_common import *
from pandas import factorize
SIZE = 1000000
indices = np.random.randint(100, size=SIZE)
"""


# --- Integer array factorization
setup = common_setup + """
int_values_uniq = np.arange(SIZE) * 100
"""
factorize_int_uniq = Benchmark("factorize(int_values_uniq)", setup,
                               start_date=START_DATE)
setup = common_setup + """
int_values_dup = (np.arange(SIZE) * 100).take(indices)
"""
factorize_int_dup = Benchmark("factorize(int_values_dup)", setup,
                              start_date=START_DATE)


# --- String array factorization
setup = common_setup + """
str_values_uniq = tm.makeStringIndex(SIZE)
"""
factorize_str_uniq = Benchmark("factorize(str_values_uniq)", setup=setup,
                               start_date=START_DATE)
setup = common_setup + """
str_values_dup = tm.makeStringIndex(SIZE).take(indices)
"""
factorize_str_dup = Benchmark("factorize(str_values_dup)", setup=setup,
                              start_date=START_DATE)
setup = common_setup + """
shortstr_4_dup = Index(np.take(['AA', 'BB', 'CC', 'DD'],
                       np.random.randint(4, size=SIZE)))
"""
factorize_shortstr_4_dup = Benchmark("factorize(shortstr_values_dup)",
                                     setup=setup, start_date=START_DATE)
setup = common_setup + """
shortstr_many_dup = tm.rands_array(2, SIZE)
"""
factorize_shortstr_many_dup = Benchmark("factorize(shortstr_many_dup)",
                                        setup=setup, start_date=START_DATE)


# --- Float array factorization
setup = common_setup + """
float_values_uniq = np.linspace(0., 1., num=SIZE) * 100
"""
factorize_float_uniq = Benchmark("factorize(float_values_uniq)", setup=setup,
                                 start_date=START_DATE)
setup = common_setup + """
float_values_dup = (np.linspace(0., 1., num=SIZE) * 100).take(indices)
"""
factorize_float_dup = Benchmark("factorize(float_values_dup)", setup,
                                start_date=START_DATE)
