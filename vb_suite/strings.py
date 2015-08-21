from vbench.api import Benchmark

common_setup = """from .pandas_vb_common import *
"""

setup = common_setup + """
import string
import itertools as IT

def make_series(letters, strlen, size):
    return Series(
        [str(x) for x in np.fromiter(IT.cycle(letters), count=size*strlen, dtype='|S1')
        .view('|S{}'.format(strlen))])

many = make_series('matchthis'+string.ascii_uppercase, strlen=19, size=10000) # 31% matches
few = make_series('matchthis'+string.ascii_uppercase*42, strlen=19, size=10000) # 1% matches
"""

strings_cat = Benchmark("many.str.cat(sep=',')", setup)
strings_title = Benchmark("many.str.title()", setup)
strings_count = Benchmark("many.str.count('matchthis')", setup)
strings_contains_many = Benchmark("many.str.contains('matchthis')", setup)
strings_contains_few = Benchmark("few.str.contains('matchthis')", setup)
strings_contains_many_noregex = Benchmark(
    "many.str.contains('matchthis', regex=False)", setup)
strings_contains_few_noregex = Benchmark(
    "few.str.contains('matchthis', regex=False)", setup)
strings_startswith = Benchmark("many.str.startswith('matchthis')", setup)
strings_endswith = Benchmark("many.str.endswith('matchthis')", setup)
strings_lower = Benchmark("many.str.lower()", setup)
strings_upper = Benchmark("many.str.upper()", setup)
strings_replace = Benchmark("many.str.replace(r'(matchthis)', r'\1\1')", setup)
strings_repeat = Benchmark(
    "many.str.repeat(list(IT.islice(IT.cycle(range(1,4)),len(many))))", setup)
strings_match = Benchmark("many.str.match(r'mat..this')", setup)
strings_extract = Benchmark("many.str.extract(r'(\w*)matchthis(\w*)')", setup)
strings_join_split = Benchmark("many.str.join(r'--').str.split('--')", setup)
strings_join_split_expand = Benchmark("many.str.join(r'--').str.split('--',expand=True)", setup)
strings_len = Benchmark("many.str.len()", setup)
strings_findall = Benchmark("many.str.findall(r'[A-Z]+')", setup)
strings_pad = Benchmark("many.str.pad(100, side='both')", setup)
strings_center = Benchmark("many.str.center(100)", setup)
strings_slice = Benchmark("many.str.slice(5,15,2)", setup)
strings_strip = Benchmark("many.str.strip('matchthis')", setup)
strings_lstrip = Benchmark("many.str.lstrip('matchthis')", setup)
strings_rstrip = Benchmark("many.str.rstrip('matchthis')", setup)
strings_get = Benchmark("many.str.get(0)", setup)

setup = setup + """
s = make_series(string.ascii_uppercase, strlen=10, size=10000).str.join('|')
"""
strings_get_dummies = Benchmark("s.str.get_dummies('|')", setup)

setup = common_setup + """
import pandas.util.testing as testing
ser = Series(testing.makeUnicodeIndex())
"""

strings_encode_decode = Benchmark("ser.str.encode('utf-8').str.decode('utf-8')", setup)
