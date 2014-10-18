from vbench.api import Benchmark
from datetime import datetime

START_DATE = datetime(2014, 10, 13)

# GH 8524

setup = """from pandas_vb_common import *
from pandas import factorize
SIZE = 1000000

int_values_uniq = np.arange(SIZE) * 100
str_values_uniq = tm.makeStringIndex(SIZE)
float_values_uniq = np.linspace(0., 1., num=SIZE) * 100

indices = np.random.randint(100, size=SIZE)
int_values_dup = int_values_uniq.take(indices)
str_values_dup = str_values_uniq.take(indices)
shortstr_values_dup = Index(np.take(['AA', 'BB', 'CC', 'DD'],
                                    np.random.randint(4, size=SIZE)))
float_values_dup = float_values_uniq.take(indices)
"""


factorize_int_uniq = Benchmark("factorize(int_values_uniq)", setup,
                               start_date=START_DATE)
factorize_int_dup = Benchmark("factorize(int_values_dup)", setup,
                              start_date=START_DATE)

factorize_str_uniq = Benchmark("factorize(str_values_uniq)", setup,
                               start_date=START_DATE)
factorize_str_dup = Benchmark("factorize(str_values_dup)", setup,
                              start_date=START_DATE)
factorize_shortstr_dup = Benchmark("factorize(shortstr_values_dup)", setup,
                                   start_date=START_DATE)

factorize_float_uniq = Benchmark("factorize(float_values_uniq)", setup,
                                 start_date=START_DATE)
factorize_float_dup = Benchmark("factorize(float_values_dup)", setup,
                                start_date=START_DATE)
