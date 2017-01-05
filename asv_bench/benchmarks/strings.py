from .pandas_vb_common import *
import string
import itertools as IT
import pandas.util.testing as testing


class StringMethods(object):
    goal_time = 0.2

    def make_series(self, letters, strlen, size):
        return Series([str(x) for x in np.fromiter(IT.cycle(letters), count=(size * strlen), dtype='|S1').view('|S{}'.format(strlen))])

    def setup(self):
        self.many = self.make_series(('matchthis' + string.ascii_uppercase), strlen=19, size=10000)
        self.few = self.make_series(('matchthis' + (string.ascii_uppercase * 42)), strlen=19, size=10000)
        self.s = self.make_series(string.ascii_uppercase, strlen=10, size=10000).str.join('|')

    def time_cat(self):
        self.many.str.cat(sep=',')

    def time_center(self):
        self.many.str.center(100)

    def time_contains_few(self):
        self.few.str.contains('matchthis')

    def time_contains_few_noregex(self):
        self.few.str.contains('matchthis', regex=False)

    def time_contains_many(self):
        self.many.str.contains('matchthis')

    def time_contains_many_noregex(self):
        self.many.str.contains('matchthis', regex=False)

    def time_count(self):
        self.many.str.count('matchthis')

    def time_endswith(self):
        self.many.str.endswith('matchthis')

    def time_extract(self):
        self.many.str.extract('(\\w*)matchthis(\\w*)')

    def time_findall(self):
        self.many.str.findall('[A-Z]+')

    def time_get(self):
        self.many.str.get(0)

    def time_join_split(self):
        self.many.str.join('--').str.split('--')

    def time_join_split_expand(self):
        self.many.str.join('--').str.split('--', expand=True)

    def time_len(self):
        self.many.str.len()

    def time_match(self):
        self.many.str.match('mat..this')

    def time_pad(self):
        self.many.str.pad(100, side='both')

    def time_repeat(self):
        self.many.str.repeat(list(IT.islice(IT.cycle(range(1, 4)), len(self.many))))

    def time_replace(self):
        self.many.str.replace('(matchthis)', '\x01\x01')

    def time_slice(self):
        self.many.str.slice(5, 15, 2)

    def time_startswith(self):
        self.many.str.startswith('matchthis')

    def time_strip(self):
        self.many.str.strip('matchthis')

    def time_rstrip(self):
        self.many.str.rstrip('matchthis')

    def time_lstrip(self):
        self.many.str.lstrip('matchthis')

    def time_title(self):
        self.many.str.title()

    def time_upper(self):
        self.many.str.upper()

    def time_lower(self):
        self.many.str.lower()

    def time_get_dummies(self):
        self.s.str.get_dummies('|')


class StringEncode(object):
    goal_time = 0.2

    def setup(self):
        self.ser = Series(testing.makeUnicodeIndex())

    def time_encode_decode(self):
        self.ser.str.encode('utf-8').str.decode('utf-8')
