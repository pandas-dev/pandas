from vbench.benchmark import Benchmark
from datetime import datetime

common_setup = """from pandas_vb_common import *
"""

#----------------------------------------------------------------------
# Panel.from_dict homogenization time

START_DATE = datetime(2011, 6, 1)

setup_same_index = common_setup + """
# create 1000 dataframes with the same index
dr = DateRange(datetime(1990,1,1), datetime(2012,1,1))
data_frames = {}
for x in xrange(1000):
   df = DataFrame({"a": [0]*len(dr), "b": [1]*len(dr),
                   "c": [2]*len(dr)}, index=dr)
   data_frames[x] = df
"""

panel_from_dict_same_index = \
    Benchmark("Panel.from_dict(data_frames)",
              setup_same_index, name='panel_from_dict_same_index',
              start_date=START_DATE, repeat=1)

setup_equiv_indexes = common_setup + """
data_frames = {}
for x in xrange(1000):
   dr = DateRange(datetime(1990,1,1), datetime(2012,1,1))
   df = DataFrame({"a": [0]*len(dr), "b": [1]*len(dr),
                   "c": [2]*len(dr)}, index=dr)
   data_frames[x] = df
"""

panel_from_dict_equiv_indexes = \
    Benchmark("Panel.from_dict(data_frames)",
              setup_equiv_indexes, name='panel_from_dict_equiv_indexes',
              start_date=START_DATE, repeat=1)

setup_all_different_indexes = common_setup + """
data_frames = {}
start = datetime(1990,1,1)
end = datetime(2012,1,1)
for x in xrange(1000):
   end += timedelta(days=1)
   dr = DateRange(start, end)
   df = DataFrame({"a": [0]*len(dr), "b": [1]*len(dr),
                   "c": [2]*len(dr)}, index=dr)
   data_frames[x] = df
"""
panel_from_dict_all_different_indexes = \
    Benchmark("Panel.from_dict(data_frames)",
              setup_all_different_indexes,
              name='panel_from_dict_all_different_indexes',
              start_date=START_DATE, repeat=1)

setup_two_different_indexes = common_setup + """
data_frames = {}
start = datetime(1990,1,1)
end = datetime(2012,1,1)
for x in xrange(1000):
   if x == 500:
       end += timedelta(days=1)
   dr = DateRange(start, end)
   df = DataFrame({"a": [0]*len(dr), "b": [1]*len(dr),
                   "c": [2]*len(dr)}, index=dr)
   data_frames[x] = df
"""
panel_from_dict_two_different_indexes = \
    Benchmark("Panel.from_dict(data_frames)",
              setup_two_different_indexes,
              name='panel_from_dict_two_different_indexes',
              start_date=START_DATE, repeat=1)
