from cStringIO import StringIO
from pandas_vb_common import *
import os
from pandas import read_csv, read_table


class read_csv_comment2(object):
    goal_time = 0.2

    def setup(self):
        self.data = ['A,B,C']
        self.data = (self.data + (['1,2,3 # comment'] * 100000))
        self.data = '\n'.join(self.data)

    def time_read_csv_comment2(self):
        read_csv(StringIO(self.data), comment='#')


class read_csv_default_converter(object):
    goal_time = 0.2

    def setup(self):
        self.data = '0.1213700904466425978256438611,0.0525708283766902484401839501,0.4174092731488769913994474336\n        0.4096341697147408700274695547,0.1587830198973579909349496119,0.1292545832485494372576795285\n        0.8323255650024565799327547210,0.9694902427379478160318626578,0.6295047811546814475747169126\n        0.4679375305798131323697930383,0.2963942381834381301075609371,0.5268936082160610157032465394\n        0.6685382761849776311890991564,0.6721207066140679753374342908,0.6519975277021627935170045020\n        '
        self.data = (self.data * 200)

    def time_read_csv_default_converter(self):
        read_csv(StringIO(self.data), sep=',', header=None, float_precision=None)


class read_csv_precise_converter(object):
    goal_time = 0.2

    def setup(self):
        self.data = '0.1213700904466425978256438611,0.0525708283766902484401839501,0.4174092731488769913994474336\n        0.4096341697147408700274695547,0.1587830198973579909349496119,0.1292545832485494372576795285\n        0.8323255650024565799327547210,0.9694902427379478160318626578,0.6295047811546814475747169126\n        0.4679375305798131323697930383,0.2963942381834381301075609371,0.5268936082160610157032465394\n        0.6685382761849776311890991564,0.6721207066140679753374342908,0.6519975277021627935170045020\n        '
        self.data = (self.data * 200)

    def time_read_csv_precise_converter(self):
        read_csv(StringIO(self.data), sep=',', header=None, float_precision='high')


class read_csv_roundtrip_converter(object):
    goal_time = 0.2

    def setup(self):
        self.data = '0.1213700904466425978256438611,0.0525708283766902484401839501,0.4174092731488769913994474336\n        0.4096341697147408700274695547,0.1587830198973579909349496119,0.1292545832485494372576795285\n        0.8323255650024565799327547210,0.9694902427379478160318626578,0.6295047811546814475747169126\n        0.4679375305798131323697930383,0.2963942381834381301075609371,0.5268936082160610157032465394\n        0.6685382761849776311890991564,0.6721207066140679753374342908,0.6519975277021627935170045020\n        '
        self.data = (self.data * 200)

    def time_read_csv_roundtrip_converter(self):
        read_csv(StringIO(self.data), sep=',', header=None, float_precision='round_trip')


class read_csv_thou_vb(object):
    goal_time = 0.2

    def setup(self):
        self.N = 10000
        self.K = 8
        self.format = (lambda x: '{:,}'.format(x))
        self.df = DataFrame((np.random.randn(self.N, self.K) * np.random.randint(100, 10000, (self.N, self.K))))
        self.df = self.df.applymap(self.format)
        self.df.to_csv('test.csv', sep='|')

    def time_read_csv_thou_vb(self):
        read_csv('test.csv', sep='|', thousands=',')

    def teardown(self):
        os.remove('test.csv')


class read_csv_vb(object):
    goal_time = 0.2

    def setup(self):
        self.N = 10000
        self.K = 8
        self.df = DataFrame((np.random.randn(self.N, self.K) * np.random.randint(100, 10000, (self.N, self.K))))
        self.df.to_csv('test.csv', sep='|')

    def time_read_csv_vb(self):
        read_csv('test.csv', sep='|')

    def teardown(self):
        os.remove('test.csv')


class read_table_multiple_date(object):
    goal_time = 0.2

    def setup(self):
        self.N = 10000
        self.K = 8
        self.data = 'KORD,19990127, 19:00:00, 18:56:00, 0.8100, 2.8100, 7.2000, 0.0000, 280.0000\n        KORD,19990127, 20:00:00, 19:56:00, 0.0100, 2.2100, 7.2000, 0.0000, 260.0000\n        KORD,19990127, 21:00:00, 20:56:00, -0.5900, 2.2100, 5.7000, 0.0000, 280.0000\n        KORD,19990127, 21:00:00, 21:18:00, -0.9900, 2.0100, 3.6000, 0.0000, 270.0000\n        KORD,19990127, 22:00:00, 21:56:00, -0.5900, 1.7100, 5.1000, 0.0000, 290.0000\n        '
        self.data = (self.data * 200)

    def time_read_table_multiple_date(self):
        read_table(StringIO(self.data), sep=',', header=None, parse_dates=[[1, 2], [1, 3]])


class read_table_multiple_date_baseline(object):
    goal_time = 0.2

    def setup(self):
        self.N = 10000
        self.K = 8
        self.data = 'KORD,19990127 19:00:00, 18:56:00, 0.8100, 2.8100, 7.2000, 0.0000, 280.0000\n        KORD,19990127 20:00:00, 19:56:00, 0.0100, 2.2100, 7.2000, 0.0000, 260.0000\n        KORD,19990127 21:00:00, 20:56:00, -0.5900, 2.2100, 5.7000, 0.0000, 280.0000\n        KORD,19990127 21:00:00, 21:18:00, -0.9900, 2.0100, 3.6000, 0.0000, 270.0000\n        KORD,19990127 22:00:00, 21:56:00, -0.5900, 1.7100, 5.1000, 0.0000, 290.0000\n        '
        self.data = (self.data * 200)

    def time_read_table_multiple_date_baseline(self):
        read_table(StringIO(self.data), sep=',', header=None, parse_dates=[1])