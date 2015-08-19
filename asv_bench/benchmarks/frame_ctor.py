from pandas_vb_common import *
try:
    from pandas.tseries.offsets import *
except:
    from pandas.core.datetools import *


class frame_ctor_dtindex_BDayx1(object):
    goal_time = 0.2

    def setup(self):

        def get_period_count(start_date, off):
            self.ten_offsets_in_days = ((start_date + (off * 10)) - start_date).days
            if (self.ten_offsets_in_days == 0):
                return 1000
            else:
                return min((9 * ((Timestamp.max - start_date).days // self.ten_offsets_in_days)), 1000)

        def get_index_for_offset(off):
            self.start_date = Timestamp('1/1/1900')
            return date_range(self.start_date, periods=min(1000, get_period_count(self.start_date, off)), freq=off)
        self.idx = get_index_for_offset(BDay(1, **{}))
        self.df = DataFrame(np.random.randn(len(self.idx), 10), index=self.idx)
        self.d = dict([(col, self.df[col]) for col in self.df.columns])

    def time_frame_ctor_dtindex_BDayx1(self):
        DataFrame(self.d)


class frame_ctor_dtindex_BDayx2(object):
    goal_time = 0.2

    def setup(self):

        def get_period_count(start_date, off):
            self.ten_offsets_in_days = ((start_date + (off * 10)) - start_date).days
            if (self.ten_offsets_in_days == 0):
                return 1000
            else:
                return min((9 * ((Timestamp.max - start_date).days // self.ten_offsets_in_days)), 1000)

        def get_index_for_offset(off):
            self.start_date = Timestamp('1/1/1900')
            return date_range(self.start_date, periods=min(1000, get_period_count(self.start_date, off)), freq=off)
        self.idx = get_index_for_offset(BDay(2, **{}))
        self.df = DataFrame(np.random.randn(len(self.idx), 10), index=self.idx)
        self.d = dict([(col, self.df[col]) for col in self.df.columns])

    def time_frame_ctor_dtindex_BDayx2(self):
        DataFrame(self.d)


class frame_ctor_dtindex_BMonthBeginx1(object):
    goal_time = 0.2

    def setup(self):

        def get_period_count(start_date, off):
            self.ten_offsets_in_days = ((start_date + (off * 10)) - start_date).days
            if (self.ten_offsets_in_days == 0):
                return 1000
            else:
                return min((9 * ((Timestamp.max - start_date).days // self.ten_offsets_in_days)), 1000)

        def get_index_for_offset(off):
            self.start_date = Timestamp('1/1/1900')
            return date_range(self.start_date, periods=min(1000, get_period_count(self.start_date, off)), freq=off)
        self.idx = get_index_for_offset(BMonthBegin(1, **{}))
        self.df = DataFrame(np.random.randn(len(self.idx), 10), index=self.idx)
        self.d = dict([(col, self.df[col]) for col in self.df.columns])

    def time_frame_ctor_dtindex_BMonthBeginx1(self):
        DataFrame(self.d)


class frame_ctor_dtindex_BMonthBeginx2(object):
    goal_time = 0.2

    def setup(self):

        def get_period_count(start_date, off):
            self.ten_offsets_in_days = ((start_date + (off * 10)) - start_date).days
            if (self.ten_offsets_in_days == 0):
                return 1000
            else:
                return min((9 * ((Timestamp.max - start_date).days // self.ten_offsets_in_days)), 1000)

        def get_index_for_offset(off):
            self.start_date = Timestamp('1/1/1900')
            return date_range(self.start_date, periods=min(1000, get_period_count(self.start_date, off)), freq=off)
        self.idx = get_index_for_offset(BMonthBegin(2, **{}))
        self.df = DataFrame(np.random.randn(len(self.idx), 10), index=self.idx)
        self.d = dict([(col, self.df[col]) for col in self.df.columns])

    def time_frame_ctor_dtindex_BMonthBeginx2(self):
        DataFrame(self.d)


class frame_ctor_dtindex_BMonthEndx1(object):
    goal_time = 0.2

    def setup(self):

        def get_period_count(start_date, off):
            self.ten_offsets_in_days = ((start_date + (off * 10)) - start_date).days
            if (self.ten_offsets_in_days == 0):
                return 1000
            else:
                return min((9 * ((Timestamp.max - start_date).days // self.ten_offsets_in_days)), 1000)

        def get_index_for_offset(off):
            self.start_date = Timestamp('1/1/1900')
            return date_range(self.start_date, periods=min(1000, get_period_count(self.start_date, off)), freq=off)
        self.idx = get_index_for_offset(BMonthEnd(1, **{}))
        self.df = DataFrame(np.random.randn(len(self.idx), 10), index=self.idx)
        self.d = dict([(col, self.df[col]) for col in self.df.columns])

    def time_frame_ctor_dtindex_BMonthEndx1(self):
        DataFrame(self.d)


class frame_ctor_dtindex_BMonthEndx2(object):
    goal_time = 0.2

    def setup(self):

        def get_period_count(start_date, off):
            self.ten_offsets_in_days = ((start_date + (off * 10)) - start_date).days
            if (self.ten_offsets_in_days == 0):
                return 1000
            else:
                return min((9 * ((Timestamp.max - start_date).days // self.ten_offsets_in_days)), 1000)

        def get_index_for_offset(off):
            self.start_date = Timestamp('1/1/1900')
            return date_range(self.start_date, periods=min(1000, get_period_count(self.start_date, off)), freq=off)
        self.idx = get_index_for_offset(BMonthEnd(2, **{}))
        self.df = DataFrame(np.random.randn(len(self.idx), 10), index=self.idx)
        self.d = dict([(col, self.df[col]) for col in self.df.columns])

    def time_frame_ctor_dtindex_BMonthEndx2(self):
        DataFrame(self.d)


class frame_ctor_dtindex_BQuarterBeginx1(object):
    goal_time = 0.2

    def setup(self):

        def get_period_count(start_date, off):
            self.ten_offsets_in_days = ((start_date + (off * 10)) - start_date).days
            if (self.ten_offsets_in_days == 0):
                return 1000
            else:
                return min((9 * ((Timestamp.max - start_date).days // self.ten_offsets_in_days)), 1000)

        def get_index_for_offset(off):
            self.start_date = Timestamp('1/1/1900')
            return date_range(self.start_date, periods=min(1000, get_period_count(self.start_date, off)), freq=off)
        self.idx = get_index_for_offset(BQuarterBegin(1, **{}))
        self.df = DataFrame(np.random.randn(len(self.idx), 10), index=self.idx)
        self.d = dict([(col, self.df[col]) for col in self.df.columns])

    def time_frame_ctor_dtindex_BQuarterBeginx1(self):
        DataFrame(self.d)


class frame_ctor_dtindex_BQuarterBeginx2(object):
    goal_time = 0.2

    def setup(self):

        def get_period_count(start_date, off):
            self.ten_offsets_in_days = ((start_date + (off * 10)) - start_date).days
            if (self.ten_offsets_in_days == 0):
                return 1000
            else:
                return min((9 * ((Timestamp.max - start_date).days // self.ten_offsets_in_days)), 1000)

        def get_index_for_offset(off):
            self.start_date = Timestamp('1/1/1900')
            return date_range(self.start_date, periods=min(1000, get_period_count(self.start_date, off)), freq=off)
        self.idx = get_index_for_offset(BQuarterBegin(2, **{}))
        self.df = DataFrame(np.random.randn(len(self.idx), 10), index=self.idx)
        self.d = dict([(col, self.df[col]) for col in self.df.columns])

    def time_frame_ctor_dtindex_BQuarterBeginx2(self):
        DataFrame(self.d)


class frame_ctor_dtindex_BQuarterEndx1(object):
    goal_time = 0.2

    def setup(self):

        def get_period_count(start_date, off):
            self.ten_offsets_in_days = ((start_date + (off * 10)) - start_date).days
            if (self.ten_offsets_in_days == 0):
                return 1000
            else:
                return min((9 * ((Timestamp.max - start_date).days // self.ten_offsets_in_days)), 1000)

        def get_index_for_offset(off):
            self.start_date = Timestamp('1/1/1900')
            return date_range(self.start_date, periods=min(1000, get_period_count(self.start_date, off)), freq=off)
        self.idx = get_index_for_offset(BQuarterEnd(1, **{}))
        self.df = DataFrame(np.random.randn(len(self.idx), 10), index=self.idx)
        self.d = dict([(col, self.df[col]) for col in self.df.columns])

    def time_frame_ctor_dtindex_BQuarterEndx1(self):
        DataFrame(self.d)


class frame_ctor_dtindex_BQuarterEndx2(object):
    goal_time = 0.2

    def setup(self):

        def get_period_count(start_date, off):
            self.ten_offsets_in_days = ((start_date + (off * 10)) - start_date).days
            if (self.ten_offsets_in_days == 0):
                return 1000
            else:
                return min((9 * ((Timestamp.max - start_date).days // self.ten_offsets_in_days)), 1000)

        def get_index_for_offset(off):
            self.start_date = Timestamp('1/1/1900')
            return date_range(self.start_date, periods=min(1000, get_period_count(self.start_date, off)), freq=off)
        self.idx = get_index_for_offset(BQuarterEnd(2, **{}))
        self.df = DataFrame(np.random.randn(len(self.idx), 10), index=self.idx)
        self.d = dict([(col, self.df[col]) for col in self.df.columns])

    def time_frame_ctor_dtindex_BQuarterEndx2(self):
        DataFrame(self.d)


class frame_ctor_dtindex_BYearBeginx1(object):
    goal_time = 0.2

    def setup(self):

        def get_period_count(start_date, off):
            self.ten_offsets_in_days = ((start_date + (off * 10)) - start_date).days
            if (self.ten_offsets_in_days == 0):
                return 1000
            else:
                return min((9 * ((Timestamp.max - start_date).days // self.ten_offsets_in_days)), 1000)

        def get_index_for_offset(off):
            self.start_date = Timestamp('1/1/1900')
            return date_range(self.start_date, periods=min(1000, get_period_count(self.start_date, off)), freq=off)
        self.idx = get_index_for_offset(BYearBegin(1, **{}))
        self.df = DataFrame(np.random.randn(len(self.idx), 10), index=self.idx)
        self.d = dict([(col, self.df[col]) for col in self.df.columns])

    def time_frame_ctor_dtindex_BYearBeginx1(self):
        DataFrame(self.d)


class frame_ctor_dtindex_BYearBeginx2(object):
    goal_time = 0.2

    def setup(self):

        def get_period_count(start_date, off):
            self.ten_offsets_in_days = ((start_date + (off * 10)) - start_date).days
            if (self.ten_offsets_in_days == 0):
                return 1000
            else:
                return min((9 * ((Timestamp.max - start_date).days // self.ten_offsets_in_days)), 1000)

        def get_index_for_offset(off):
            self.start_date = Timestamp('1/1/1900')
            return date_range(self.start_date, periods=min(1000, get_period_count(self.start_date, off)), freq=off)
        self.idx = get_index_for_offset(BYearBegin(2, **{}))
        self.df = DataFrame(np.random.randn(len(self.idx), 10), index=self.idx)
        self.d = dict([(col, self.df[col]) for col in self.df.columns])

    def time_frame_ctor_dtindex_BYearBeginx2(self):
        DataFrame(self.d)


class frame_ctor_dtindex_BYearEndx1(object):
    goal_time = 0.2

    def setup(self):

        def get_period_count(start_date, off):
            self.ten_offsets_in_days = ((start_date + (off * 10)) - start_date).days
            if (self.ten_offsets_in_days == 0):
                return 1000
            else:
                return min((9 * ((Timestamp.max - start_date).days // self.ten_offsets_in_days)), 1000)

        def get_index_for_offset(off):
            self.start_date = Timestamp('1/1/1900')
            return date_range(self.start_date, periods=min(1000, get_period_count(self.start_date, off)), freq=off)
        self.idx = get_index_for_offset(BYearEnd(1, **{}))
        self.df = DataFrame(np.random.randn(len(self.idx), 10), index=self.idx)
        self.d = dict([(col, self.df[col]) for col in self.df.columns])

    def time_frame_ctor_dtindex_BYearEndx1(self):
        DataFrame(self.d)


class frame_ctor_dtindex_BYearEndx2(object):
    goal_time = 0.2

    def setup(self):

        def get_period_count(start_date, off):
            self.ten_offsets_in_days = ((start_date + (off * 10)) - start_date).days
            if (self.ten_offsets_in_days == 0):
                return 1000
            else:
                return min((9 * ((Timestamp.max - start_date).days // self.ten_offsets_in_days)), 1000)

        def get_index_for_offset(off):
            self.start_date = Timestamp('1/1/1900')
            return date_range(self.start_date, periods=min(1000, get_period_count(self.start_date, off)), freq=off)
        self.idx = get_index_for_offset(BYearEnd(2, **{}))
        self.df = DataFrame(np.random.randn(len(self.idx), 10), index=self.idx)
        self.d = dict([(col, self.df[col]) for col in self.df.columns])

    def time_frame_ctor_dtindex_BYearEndx2(self):
        DataFrame(self.d)


class frame_ctor_dtindex_BusinessDayx1(object):
    goal_time = 0.2

    def setup(self):

        def get_period_count(start_date, off):
            self.ten_offsets_in_days = ((start_date + (off * 10)) - start_date).days
            if (self.ten_offsets_in_days == 0):
                return 1000
            else:
                return min((9 * ((Timestamp.max - start_date).days // self.ten_offsets_in_days)), 1000)

        def get_index_for_offset(off):
            self.start_date = Timestamp('1/1/1900')
            return date_range(self.start_date, periods=min(1000, get_period_count(self.start_date, off)), freq=off)
        self.idx = get_index_for_offset(BusinessDay(1, **{}))
        self.df = DataFrame(np.random.randn(len(self.idx), 10), index=self.idx)
        self.d = dict([(col, self.df[col]) for col in self.df.columns])

    def time_frame_ctor_dtindex_BusinessDayx1(self):
        DataFrame(self.d)


class frame_ctor_dtindex_BusinessDayx2(object):
    goal_time = 0.2

    def setup(self):

        def get_period_count(start_date, off):
            self.ten_offsets_in_days = ((start_date + (off * 10)) - start_date).days
            if (self.ten_offsets_in_days == 0):
                return 1000
            else:
                return min((9 * ((Timestamp.max - start_date).days // self.ten_offsets_in_days)), 1000)

        def get_index_for_offset(off):
            self.start_date = Timestamp('1/1/1900')
            return date_range(self.start_date, periods=min(1000, get_period_count(self.start_date, off)), freq=off)
        self.idx = get_index_for_offset(BusinessDay(2, **{}))
        self.df = DataFrame(np.random.randn(len(self.idx), 10), index=self.idx)
        self.d = dict([(col, self.df[col]) for col in self.df.columns])

    def time_frame_ctor_dtindex_BusinessDayx2(self):
        DataFrame(self.d)


class frame_ctor_dtindex_BusinessHourx1(object):
    goal_time = 0.2

    def setup(self):

        def get_period_count(start_date, off):
            self.ten_offsets_in_days = ((start_date + (off * 10)) - start_date).days
            if (self.ten_offsets_in_days == 0):
                return 1000
            else:
                return min((9 * ((Timestamp.max - start_date).days // self.ten_offsets_in_days)), 1000)

        def get_index_for_offset(off):
            self.start_date = Timestamp('1/1/1900')
            return date_range(self.start_date, periods=min(1000, get_period_count(self.start_date, off)), freq=off)
        self.idx = get_index_for_offset(BusinessHour(1, **{}))
        self.df = DataFrame(np.random.randn(len(self.idx), 10), index=self.idx)
        self.d = dict([(col, self.df[col]) for col in self.df.columns])

    def time_frame_ctor_dtindex_BusinessHourx1(self):
        DataFrame(self.d)


class frame_ctor_dtindex_BusinessHourx2(object):
    goal_time = 0.2

    def setup(self):

        def get_period_count(start_date, off):
            self.ten_offsets_in_days = ((start_date + (off * 10)) - start_date).days
            if (self.ten_offsets_in_days == 0):
                return 1000
            else:
                return min((9 * ((Timestamp.max - start_date).days // self.ten_offsets_in_days)), 1000)

        def get_index_for_offset(off):
            self.start_date = Timestamp('1/1/1900')
            return date_range(self.start_date, periods=min(1000, get_period_count(self.start_date, off)), freq=off)
        self.idx = get_index_for_offset(BusinessHour(2, **{}))
        self.df = DataFrame(np.random.randn(len(self.idx), 10), index=self.idx)
        self.d = dict([(col, self.df[col]) for col in self.df.columns])

    def time_frame_ctor_dtindex_BusinessHourx2(self):
        DataFrame(self.d)


class frame_ctor_dtindex_CBMonthBeginx1(object):
    goal_time = 0.2

    def setup(self):

        def get_period_count(start_date, off):
            self.ten_offsets_in_days = ((start_date + (off * 10)) - start_date).days
            if (self.ten_offsets_in_days == 0):
                return 1000
            else:
                return min((9 * ((Timestamp.max - start_date).days // self.ten_offsets_in_days)), 1000)

        def get_index_for_offset(off):
            self.start_date = Timestamp('1/1/1900')
            return date_range(self.start_date, periods=min(1000, get_period_count(self.start_date, off)), freq=off)
        self.idx = get_index_for_offset(CBMonthBegin(1, **{}))
        self.df = DataFrame(np.random.randn(len(self.idx), 10), index=self.idx)
        self.d = dict([(col, self.df[col]) for col in self.df.columns])

    def time_frame_ctor_dtindex_CBMonthBeginx1(self):
        DataFrame(self.d)


class frame_ctor_dtindex_CBMonthBeginx2(object):
    goal_time = 0.2

    def setup(self):

        def get_period_count(start_date, off):
            self.ten_offsets_in_days = ((start_date + (off * 10)) - start_date).days
            if (self.ten_offsets_in_days == 0):
                return 1000
            else:
                return min((9 * ((Timestamp.max - start_date).days // self.ten_offsets_in_days)), 1000)

        def get_index_for_offset(off):
            self.start_date = Timestamp('1/1/1900')
            return date_range(self.start_date, periods=min(1000, get_period_count(self.start_date, off)), freq=off)
        self.idx = get_index_for_offset(CBMonthBegin(2, **{}))
        self.df = DataFrame(np.random.randn(len(self.idx), 10), index=self.idx)
        self.d = dict([(col, self.df[col]) for col in self.df.columns])

    def time_frame_ctor_dtindex_CBMonthBeginx2(self):
        DataFrame(self.d)


class frame_ctor_dtindex_CBMonthEndx1(object):
    goal_time = 0.2

    def setup(self):

        def get_period_count(start_date, off):
            self.ten_offsets_in_days = ((start_date + (off * 10)) - start_date).days
            if (self.ten_offsets_in_days == 0):
                return 1000
            else:
                return min((9 * ((Timestamp.max - start_date).days // self.ten_offsets_in_days)), 1000)

        def get_index_for_offset(off):
            self.start_date = Timestamp('1/1/1900')
            return date_range(self.start_date, periods=min(1000, get_period_count(self.start_date, off)), freq=off)
        self.idx = get_index_for_offset(CBMonthEnd(1, **{}))
        self.df = DataFrame(np.random.randn(len(self.idx), 10), index=self.idx)
        self.d = dict([(col, self.df[col]) for col in self.df.columns])

    def time_frame_ctor_dtindex_CBMonthEndx1(self):
        DataFrame(self.d)


class frame_ctor_dtindex_CBMonthEndx2(object):
    goal_time = 0.2

    def setup(self):

        def get_period_count(start_date, off):
            self.ten_offsets_in_days = ((start_date + (off * 10)) - start_date).days
            if (self.ten_offsets_in_days == 0):
                return 1000
            else:
                return min((9 * ((Timestamp.max - start_date).days // self.ten_offsets_in_days)), 1000)

        def get_index_for_offset(off):
            self.start_date = Timestamp('1/1/1900')
            return date_range(self.start_date, periods=min(1000, get_period_count(self.start_date, off)), freq=off)
        self.idx = get_index_for_offset(CBMonthEnd(2, **{}))
        self.df = DataFrame(np.random.randn(len(self.idx), 10), index=self.idx)
        self.d = dict([(col, self.df[col]) for col in self.df.columns])

    def time_frame_ctor_dtindex_CBMonthEndx2(self):
        DataFrame(self.d)


class frame_ctor_dtindex_CDayx1(object):
    goal_time = 0.2

    def setup(self):

        def get_period_count(start_date, off):
            self.ten_offsets_in_days = ((start_date + (off * 10)) - start_date).days
            if (self.ten_offsets_in_days == 0):
                return 1000
            else:
                return min((9 * ((Timestamp.max - start_date).days // self.ten_offsets_in_days)), 1000)

        def get_index_for_offset(off):
            self.start_date = Timestamp('1/1/1900')
            return date_range(self.start_date, periods=min(1000, get_period_count(self.start_date, off)), freq=off)
        self.idx = get_index_for_offset(CDay(1, **{}))
        self.df = DataFrame(np.random.randn(len(self.idx), 10), index=self.idx)
        self.d = dict([(col, self.df[col]) for col in self.df.columns])

    def time_frame_ctor_dtindex_CDayx1(self):
        DataFrame(self.d)


class frame_ctor_dtindex_CDayx2(object):
    goal_time = 0.2

    def setup(self):

        def get_period_count(start_date, off):
            self.ten_offsets_in_days = ((start_date + (off * 10)) - start_date).days
            if (self.ten_offsets_in_days == 0):
                return 1000
            else:
                return min((9 * ((Timestamp.max - start_date).days // self.ten_offsets_in_days)), 1000)

        def get_index_for_offset(off):
            self.start_date = Timestamp('1/1/1900')
            return date_range(self.start_date, periods=min(1000, get_period_count(self.start_date, off)), freq=off)
        self.idx = get_index_for_offset(CDay(2, **{}))
        self.df = DataFrame(np.random.randn(len(self.idx), 10), index=self.idx)
        self.d = dict([(col, self.df[col]) for col in self.df.columns])

    def time_frame_ctor_dtindex_CDayx2(self):
        DataFrame(self.d)


class frame_ctor_dtindex_CustomBusinessDayx1(object):
    goal_time = 0.2

    def setup(self):

        def get_period_count(start_date, off):
            self.ten_offsets_in_days = ((start_date + (off * 10)) - start_date).days
            if (self.ten_offsets_in_days == 0):
                return 1000
            else:
                return min((9 * ((Timestamp.max - start_date).days // self.ten_offsets_in_days)), 1000)

        def get_index_for_offset(off):
            self.start_date = Timestamp('1/1/1900')
            return date_range(self.start_date, periods=min(1000, get_period_count(self.start_date, off)), freq=off)
        self.idx = get_index_for_offset(CustomBusinessDay(1, **{}))
        self.df = DataFrame(np.random.randn(len(self.idx), 10), index=self.idx)
        self.d = dict([(col, self.df[col]) for col in self.df.columns])

    def time_frame_ctor_dtindex_CustomBusinessDayx1(self):
        DataFrame(self.d)


class frame_ctor_dtindex_CustomBusinessDayx2(object):
    goal_time = 0.2

    def setup(self):

        def get_period_count(start_date, off):
            self.ten_offsets_in_days = ((start_date + (off * 10)) - start_date).days
            if (self.ten_offsets_in_days == 0):
                return 1000
            else:
                return min((9 * ((Timestamp.max - start_date).days // self.ten_offsets_in_days)), 1000)

        def get_index_for_offset(off):
            self.start_date = Timestamp('1/1/1900')
            return date_range(self.start_date, periods=min(1000, get_period_count(self.start_date, off)), freq=off)
        self.idx = get_index_for_offset(CustomBusinessDay(2, **{}))
        self.df = DataFrame(np.random.randn(len(self.idx), 10), index=self.idx)
        self.d = dict([(col, self.df[col]) for col in self.df.columns])

    def time_frame_ctor_dtindex_CustomBusinessDayx2(self):
        DataFrame(self.d)


class frame_ctor_dtindex_DateOffsetx1(object):
    goal_time = 0.2

    def setup(self):

        def get_period_count(start_date, off):
            self.ten_offsets_in_days = ((start_date + (off * 10)) - start_date).days
            if (self.ten_offsets_in_days == 0):
                return 1000
            else:
                return min((9 * ((Timestamp.max - start_date).days // self.ten_offsets_in_days)), 1000)

        def get_index_for_offset(off):
            self.start_date = Timestamp('1/1/1900')
            return date_range(self.start_date, periods=min(1000, get_period_count(self.start_date, off)), freq=off)
        self.idx = get_index_for_offset(DateOffset(1, **{}))
        self.df = DataFrame(np.random.randn(len(self.idx), 10), index=self.idx)
        self.d = dict([(col, self.df[col]) for col in self.df.columns])

    def time_frame_ctor_dtindex_DateOffsetx1(self):
        DataFrame(self.d)


class frame_ctor_dtindex_DateOffsetx2(object):
    goal_time = 0.2

    def setup(self):

        def get_period_count(start_date, off):
            self.ten_offsets_in_days = ((start_date + (off * 10)) - start_date).days
            if (self.ten_offsets_in_days == 0):
                return 1000
            else:
                return min((9 * ((Timestamp.max - start_date).days // self.ten_offsets_in_days)), 1000)

        def get_index_for_offset(off):
            self.start_date = Timestamp('1/1/1900')
            return date_range(self.start_date, periods=min(1000, get_period_count(self.start_date, off)), freq=off)
        self.idx = get_index_for_offset(DateOffset(2, **{}))
        self.df = DataFrame(np.random.randn(len(self.idx), 10), index=self.idx)
        self.d = dict([(col, self.df[col]) for col in self.df.columns])

    def time_frame_ctor_dtindex_DateOffsetx2(self):
        DataFrame(self.d)


class frame_ctor_dtindex_Dayx1(object):
    goal_time = 0.2

    def setup(self):

        def get_period_count(start_date, off):
            self.ten_offsets_in_days = ((start_date + (off * 10)) - start_date).days
            if (self.ten_offsets_in_days == 0):
                return 1000
            else:
                return min((9 * ((Timestamp.max - start_date).days // self.ten_offsets_in_days)), 1000)

        def get_index_for_offset(off):
            self.start_date = Timestamp('1/1/1900')
            return date_range(self.start_date, periods=min(1000, get_period_count(self.start_date, off)), freq=off)
        self.idx = get_index_for_offset(Day(1, **{}))
        self.df = DataFrame(np.random.randn(len(self.idx), 10), index=self.idx)
        self.d = dict([(col, self.df[col]) for col in self.df.columns])

    def time_frame_ctor_dtindex_Dayx1(self):
        DataFrame(self.d)


class frame_ctor_dtindex_Dayx2(object):
    goal_time = 0.2

    def setup(self):

        def get_period_count(start_date, off):
            self.ten_offsets_in_days = ((start_date + (off * 10)) - start_date).days
            if (self.ten_offsets_in_days == 0):
                return 1000
            else:
                return min((9 * ((Timestamp.max - start_date).days // self.ten_offsets_in_days)), 1000)

        def get_index_for_offset(off):
            self.start_date = Timestamp('1/1/1900')
            return date_range(self.start_date, periods=min(1000, get_period_count(self.start_date, off)), freq=off)
        self.idx = get_index_for_offset(Day(2, **{}))
        self.df = DataFrame(np.random.randn(len(self.idx), 10), index=self.idx)
        self.d = dict([(col, self.df[col]) for col in self.df.columns])

    def time_frame_ctor_dtindex_Dayx2(self):
        DataFrame(self.d)


class frame_ctor_dtindex_Easterx1(object):
    goal_time = 0.2

    def setup(self):

        def get_period_count(start_date, off):
            self.ten_offsets_in_days = ((start_date + (off * 10)) - start_date).days
            if (self.ten_offsets_in_days == 0):
                return 1000
            else:
                return min((9 * ((Timestamp.max - start_date).days // self.ten_offsets_in_days)), 1000)

        def get_index_for_offset(off):
            self.start_date = Timestamp('1/1/1900')
            return date_range(self.start_date, periods=min(1000, get_period_count(self.start_date, off)), freq=off)
        self.idx = get_index_for_offset(Easter(1, **{}))
        self.df = DataFrame(np.random.randn(len(self.idx), 10), index=self.idx)
        self.d = dict([(col, self.df[col]) for col in self.df.columns])

    def time_frame_ctor_dtindex_Easterx1(self):
        DataFrame(self.d)


class frame_ctor_dtindex_Easterx2(object):
    goal_time = 0.2

    def setup(self):

        def get_period_count(start_date, off):
            self.ten_offsets_in_days = ((start_date + (off * 10)) - start_date).days
            if (self.ten_offsets_in_days == 0):
                return 1000
            else:
                return min((9 * ((Timestamp.max - start_date).days // self.ten_offsets_in_days)), 1000)

        def get_index_for_offset(off):
            self.start_date = Timestamp('1/1/1900')
            return date_range(self.start_date, periods=min(1000, get_period_count(self.start_date, off)), freq=off)
        self.idx = get_index_for_offset(Easter(2, **{}))
        self.df = DataFrame(np.random.randn(len(self.idx), 10), index=self.idx)
        self.d = dict([(col, self.df[col]) for col in self.df.columns])

    def time_frame_ctor_dtindex_Easterx2(self):
        DataFrame(self.d)


class frame_ctor_dtindex_FY5253Quarterx1__variation_last(object):
    goal_time = 0.2

    def setup(self):

        def get_period_count(start_date, off):
            self.ten_offsets_in_days = ((start_date + (off * 10)) - start_date).days
            if (self.ten_offsets_in_days == 0):
                return 1000
            else:
                return min((9 * ((Timestamp.max - start_date).days // self.ten_offsets_in_days)), 1000)

        def get_index_for_offset(off):
            self.start_date = Timestamp('1/1/1900')
            return date_range(self.start_date, periods=min(1000, get_period_count(self.start_date, off)), freq=off)
        self.idx = get_index_for_offset(FY5253Quarter(1, **{'startingMonth': 1, 'qtr_with_extra_week': 1, 'weekday': 1, 'variation': 'last', }))
        self.df = DataFrame(np.random.randn(len(self.idx), 10), index=self.idx)
        self.d = dict([(col, self.df[col]) for col in self.df.columns])

    def time_frame_ctor_dtindex_FY5253Quarterx1__variation_last(self):
        DataFrame(self.d)


class frame_ctor_dtindex_FY5253Quarterx1__variation_nearest(object):
    goal_time = 0.2

    def setup(self):

        def get_period_count(start_date, off):
            self.ten_offsets_in_days = ((start_date + (off * 10)) - start_date).days
            if (self.ten_offsets_in_days == 0):
                return 1000
            else:
                return min((9 * ((Timestamp.max - start_date).days // self.ten_offsets_in_days)), 1000)

        def get_index_for_offset(off):
            self.start_date = Timestamp('1/1/1900')
            return date_range(self.start_date, periods=min(1000, get_period_count(self.start_date, off)), freq=off)
        self.idx = get_index_for_offset(FY5253Quarter(1, **{'startingMonth': 1, 'qtr_with_extra_week': 1, 'weekday': 1, 'variation': 'nearest', }))
        self.df = DataFrame(np.random.randn(len(self.idx), 10), index=self.idx)
        self.d = dict([(col, self.df[col]) for col in self.df.columns])

    def time_frame_ctor_dtindex_FY5253Quarterx1__variation_nearest(self):
        DataFrame(self.d)


class frame_ctor_dtindex_FY5253Quarterx2__variation_last(object):
    goal_time = 0.2

    def setup(self):

        def get_period_count(start_date, off):
            self.ten_offsets_in_days = ((start_date + (off * 10)) - start_date).days
            if (self.ten_offsets_in_days == 0):
                return 1000
            else:
                return min((9 * ((Timestamp.max - start_date).days // self.ten_offsets_in_days)), 1000)

        def get_index_for_offset(off):
            self.start_date = Timestamp('1/1/1900')
            return date_range(self.start_date, periods=min(1000, get_period_count(self.start_date, off)), freq=off)
        self.idx = get_index_for_offset(FY5253Quarter(2, **{'startingMonth': 1, 'qtr_with_extra_week': 1, 'weekday': 1, 'variation': 'last', }))
        self.df = DataFrame(np.random.randn(len(self.idx), 10), index=self.idx)
        self.d = dict([(col, self.df[col]) for col in self.df.columns])

    def time_frame_ctor_dtindex_FY5253Quarterx2__variation_last(self):
        DataFrame(self.d)


class frame_ctor_dtindex_FY5253Quarterx2__variation_nearest(object):
    goal_time = 0.2

    def setup(self):

        def get_period_count(start_date, off):
            self.ten_offsets_in_days = ((start_date + (off * 10)) - start_date).days
            if (self.ten_offsets_in_days == 0):
                return 1000
            else:
                return min((9 * ((Timestamp.max - start_date).days // self.ten_offsets_in_days)), 1000)

        def get_index_for_offset(off):
            self.start_date = Timestamp('1/1/1900')
            return date_range(self.start_date, periods=min(1000, get_period_count(self.start_date, off)), freq=off)
        self.idx = get_index_for_offset(FY5253Quarter(2, **{'startingMonth': 1, 'qtr_with_extra_week': 1, 'weekday': 1, 'variation': 'nearest', }))
        self.df = DataFrame(np.random.randn(len(self.idx), 10), index=self.idx)
        self.d = dict([(col, self.df[col]) for col in self.df.columns])

    def time_frame_ctor_dtindex_FY5253Quarterx2__variation_nearest(self):
        DataFrame(self.d)


class frame_ctor_dtindex_FY5253x1__variation_last(object):
    goal_time = 0.2

    def setup(self):

        def get_period_count(start_date, off):
            self.ten_offsets_in_days = ((start_date + (off * 10)) - start_date).days
            if (self.ten_offsets_in_days == 0):
                return 1000
            else:
                return min((9 * ((Timestamp.max - start_date).days // self.ten_offsets_in_days)), 1000)

        def get_index_for_offset(off):
            self.start_date = Timestamp('1/1/1900')
            return date_range(self.start_date, periods=min(1000, get_period_count(self.start_date, off)), freq=off)
        self.idx = get_index_for_offset(FY5253(1, **{'startingMonth': 1, 'weekday': 1, 'variation': 'last', }))
        self.df = DataFrame(np.random.randn(len(self.idx), 10), index=self.idx)
        self.d = dict([(col, self.df[col]) for col in self.df.columns])

    def time_frame_ctor_dtindex_FY5253x1__variation_last(self):
        DataFrame(self.d)


class frame_ctor_dtindex_FY5253x1__variation_nearest(object):
    goal_time = 0.2

    def setup(self):

        def get_period_count(start_date, off):
            self.ten_offsets_in_days = ((start_date + (off * 10)) - start_date).days
            if (self.ten_offsets_in_days == 0):
                return 1000
            else:
                return min((9 * ((Timestamp.max - start_date).days // self.ten_offsets_in_days)), 1000)

        def get_index_for_offset(off):
            self.start_date = Timestamp('1/1/1900')
            return date_range(self.start_date, periods=min(1000, get_period_count(self.start_date, off)), freq=off)
        self.idx = get_index_for_offset(FY5253(1, **{'startingMonth': 1, 'weekday': 1, 'variation': 'nearest', }))
        self.df = DataFrame(np.random.randn(len(self.idx), 10), index=self.idx)
        self.d = dict([(col, self.df[col]) for col in self.df.columns])

    def time_frame_ctor_dtindex_FY5253x1__variation_nearest(self):
        DataFrame(self.d)


class frame_ctor_dtindex_FY5253x2__variation_last(object):
    goal_time = 0.2

    def setup(self):

        def get_period_count(start_date, off):
            self.ten_offsets_in_days = ((start_date + (off * 10)) - start_date).days
            if (self.ten_offsets_in_days == 0):
                return 1000
            else:
                return min((9 * ((Timestamp.max - start_date).days // self.ten_offsets_in_days)), 1000)

        def get_index_for_offset(off):
            self.start_date = Timestamp('1/1/1900')
            return date_range(self.start_date, periods=min(1000, get_period_count(self.start_date, off)), freq=off)
        self.idx = get_index_for_offset(FY5253(2, **{'startingMonth': 1, 'weekday': 1, 'variation': 'last', }))
        self.df = DataFrame(np.random.randn(len(self.idx), 10), index=self.idx)
        self.d = dict([(col, self.df[col]) for col in self.df.columns])

    def time_frame_ctor_dtindex_FY5253x2__variation_last(self):
        DataFrame(self.d)


class frame_ctor_dtindex_FY5253x2__variation_nearest(object):
    goal_time = 0.2

    def setup(self):

        def get_period_count(start_date, off):
            self.ten_offsets_in_days = ((start_date + (off * 10)) - start_date).days
            if (self.ten_offsets_in_days == 0):
                return 1000
            else:
                return min((9 * ((Timestamp.max - start_date).days // self.ten_offsets_in_days)), 1000)

        def get_index_for_offset(off):
            self.start_date = Timestamp('1/1/1900')
            return date_range(self.start_date, periods=min(1000, get_period_count(self.start_date, off)), freq=off)
        self.idx = get_index_for_offset(FY5253(2, **{'startingMonth': 1, 'weekday': 1, 'variation': 'nearest', }))
        self.df = DataFrame(np.random.randn(len(self.idx), 10), index=self.idx)
        self.d = dict([(col, self.df[col]) for col in self.df.columns])

    def time_frame_ctor_dtindex_FY5253x2__variation_nearest(self):
        DataFrame(self.d)


class frame_ctor_dtindex_Hourx1(object):
    goal_time = 0.2

    def setup(self):

        def get_period_count(start_date, off):
            self.ten_offsets_in_days = ((start_date + (off * 10)) - start_date).days
            if (self.ten_offsets_in_days == 0):
                return 1000
            else:
                return min((9 * ((Timestamp.max - start_date).days // self.ten_offsets_in_days)), 1000)

        def get_index_for_offset(off):
            self.start_date = Timestamp('1/1/1900')
            return date_range(self.start_date, periods=min(1000, get_period_count(self.start_date, off)), freq=off)
        self.idx = get_index_for_offset(Hour(1, **{}))
        self.df = DataFrame(np.random.randn(len(self.idx), 10), index=self.idx)
        self.d = dict([(col, self.df[col]) for col in self.df.columns])

    def time_frame_ctor_dtindex_Hourx1(self):
        DataFrame(self.d)


class frame_ctor_dtindex_Hourx2(object):
    goal_time = 0.2

    def setup(self):

        def get_period_count(start_date, off):
            self.ten_offsets_in_days = ((start_date + (off * 10)) - start_date).days
            if (self.ten_offsets_in_days == 0):
                return 1000
            else:
                return min((9 * ((Timestamp.max - start_date).days // self.ten_offsets_in_days)), 1000)

        def get_index_for_offset(off):
            self.start_date = Timestamp('1/1/1900')
            return date_range(self.start_date, periods=min(1000, get_period_count(self.start_date, off)), freq=off)
        self.idx = get_index_for_offset(Hour(2, **{}))
        self.df = DataFrame(np.random.randn(len(self.idx), 10), index=self.idx)
        self.d = dict([(col, self.df[col]) for col in self.df.columns])

    def time_frame_ctor_dtindex_Hourx2(self):
        DataFrame(self.d)


class frame_ctor_dtindex_LastWeekOfMonthx1(object):
    goal_time = 0.2

    def setup(self):

        def get_period_count(start_date, off):
            self.ten_offsets_in_days = ((start_date + (off * 10)) - start_date).days
            if (self.ten_offsets_in_days == 0):
                return 1000
            else:
                return min((9 * ((Timestamp.max - start_date).days // self.ten_offsets_in_days)), 1000)

        def get_index_for_offset(off):
            self.start_date = Timestamp('1/1/1900')
            return date_range(self.start_date, periods=min(1000, get_period_count(self.start_date, off)), freq=off)
        self.idx = get_index_for_offset(LastWeekOfMonth(1, **{'week': 1, 'weekday': 1, }))
        self.df = DataFrame(np.random.randn(len(self.idx), 10), index=self.idx)
        self.d = dict([(col, self.df[col]) for col in self.df.columns])

    def time_frame_ctor_dtindex_LastWeekOfMonthx1(self):
        DataFrame(self.d)


class frame_ctor_dtindex_LastWeekOfMonthx2(object):
    goal_time = 0.2

    def setup(self):

        def get_period_count(start_date, off):
            self.ten_offsets_in_days = ((start_date + (off * 10)) - start_date).days
            if (self.ten_offsets_in_days == 0):
                return 1000
            else:
                return min((9 * ((Timestamp.max - start_date).days // self.ten_offsets_in_days)), 1000)

        def get_index_for_offset(off):
            self.start_date = Timestamp('1/1/1900')
            return date_range(self.start_date, periods=min(1000, get_period_count(self.start_date, off)), freq=off)
        self.idx = get_index_for_offset(LastWeekOfMonth(2, **{'week': 1, 'weekday': 1, }))
        self.df = DataFrame(np.random.randn(len(self.idx), 10), index=self.idx)
        self.d = dict([(col, self.df[col]) for col in self.df.columns])

    def time_frame_ctor_dtindex_LastWeekOfMonthx2(self):
        DataFrame(self.d)


class frame_ctor_dtindex_Microx1(object):
    goal_time = 0.2

    def setup(self):

        def get_period_count(start_date, off):
            self.ten_offsets_in_days = ((start_date + (off * 10)) - start_date).days
            if (self.ten_offsets_in_days == 0):
                return 1000
            else:
                return min((9 * ((Timestamp.max - start_date).days // self.ten_offsets_in_days)), 1000)

        def get_index_for_offset(off):
            self.start_date = Timestamp('1/1/1900')
            return date_range(self.start_date, periods=min(1000, get_period_count(self.start_date, off)), freq=off)
        self.idx = get_index_for_offset(Micro(1, **{}))
        self.df = DataFrame(np.random.randn(len(self.idx), 10), index=self.idx)
        self.d = dict([(col, self.df[col]) for col in self.df.columns])

    def time_frame_ctor_dtindex_Microx1(self):
        DataFrame(self.d)


class frame_ctor_dtindex_Microx2(object):
    goal_time = 0.2

    def setup(self):

        def get_period_count(start_date, off):
            self.ten_offsets_in_days = ((start_date + (off * 10)) - start_date).days
            if (self.ten_offsets_in_days == 0):
                return 1000
            else:
                return min((9 * ((Timestamp.max - start_date).days // self.ten_offsets_in_days)), 1000)

        def get_index_for_offset(off):
            self.start_date = Timestamp('1/1/1900')
            return date_range(self.start_date, periods=min(1000, get_period_count(self.start_date, off)), freq=off)
        self.idx = get_index_for_offset(Micro(2, **{}))
        self.df = DataFrame(np.random.randn(len(self.idx), 10), index=self.idx)
        self.d = dict([(col, self.df[col]) for col in self.df.columns])

    def time_frame_ctor_dtindex_Microx2(self):
        DataFrame(self.d)


class frame_ctor_dtindex_Millix1(object):
    goal_time = 0.2

    def setup(self):

        def get_period_count(start_date, off):
            self.ten_offsets_in_days = ((start_date + (off * 10)) - start_date).days
            if (self.ten_offsets_in_days == 0):
                return 1000
            else:
                return min((9 * ((Timestamp.max - start_date).days // self.ten_offsets_in_days)), 1000)

        def get_index_for_offset(off):
            self.start_date = Timestamp('1/1/1900')
            return date_range(self.start_date, periods=min(1000, get_period_count(self.start_date, off)), freq=off)
        self.idx = get_index_for_offset(Milli(1, **{}))
        self.df = DataFrame(np.random.randn(len(self.idx), 10), index=self.idx)
        self.d = dict([(col, self.df[col]) for col in self.df.columns])

    def time_frame_ctor_dtindex_Millix1(self):
        DataFrame(self.d)


class frame_ctor_dtindex_Millix2(object):
    goal_time = 0.2

    def setup(self):

        def get_period_count(start_date, off):
            self.ten_offsets_in_days = ((start_date + (off * 10)) - start_date).days
            if (self.ten_offsets_in_days == 0):
                return 1000
            else:
                return min((9 * ((Timestamp.max - start_date).days // self.ten_offsets_in_days)), 1000)

        def get_index_for_offset(off):
            self.start_date = Timestamp('1/1/1900')
            return date_range(self.start_date, periods=min(1000, get_period_count(self.start_date, off)), freq=off)
        self.idx = get_index_for_offset(Milli(2, **{}))
        self.df = DataFrame(np.random.randn(len(self.idx), 10), index=self.idx)
        self.d = dict([(col, self.df[col]) for col in self.df.columns])

    def time_frame_ctor_dtindex_Millix2(self):
        DataFrame(self.d)


class frame_ctor_dtindex_Minutex1(object):
    goal_time = 0.2

    def setup(self):

        def get_period_count(start_date, off):
            self.ten_offsets_in_days = ((start_date + (off * 10)) - start_date).days
            if (self.ten_offsets_in_days == 0):
                return 1000
            else:
                return min((9 * ((Timestamp.max - start_date).days // self.ten_offsets_in_days)), 1000)

        def get_index_for_offset(off):
            self.start_date = Timestamp('1/1/1900')
            return date_range(self.start_date, periods=min(1000, get_period_count(self.start_date, off)), freq=off)
        self.idx = get_index_for_offset(Minute(1, **{}))
        self.df = DataFrame(np.random.randn(len(self.idx), 10), index=self.idx)
        self.d = dict([(col, self.df[col]) for col in self.df.columns])

    def time_frame_ctor_dtindex_Minutex1(self):
        DataFrame(self.d)


class frame_ctor_dtindex_Minutex2(object):
    goal_time = 0.2

    def setup(self):

        def get_period_count(start_date, off):
            self.ten_offsets_in_days = ((start_date + (off * 10)) - start_date).days
            if (self.ten_offsets_in_days == 0):
                return 1000
            else:
                return min((9 * ((Timestamp.max - start_date).days // self.ten_offsets_in_days)), 1000)

        def get_index_for_offset(off):
            self.start_date = Timestamp('1/1/1900')
            return date_range(self.start_date, periods=min(1000, get_period_count(self.start_date, off)), freq=off)
        self.idx = get_index_for_offset(Minute(2, **{}))
        self.df = DataFrame(np.random.randn(len(self.idx), 10), index=self.idx)
        self.d = dict([(col, self.df[col]) for col in self.df.columns])

    def time_frame_ctor_dtindex_Minutex2(self):
        DataFrame(self.d)


class frame_ctor_dtindex_MonthBeginx1(object):
    goal_time = 0.2

    def setup(self):

        def get_period_count(start_date, off):
            self.ten_offsets_in_days = ((start_date + (off * 10)) - start_date).days
            if (self.ten_offsets_in_days == 0):
                return 1000
            else:
                return min((9 * ((Timestamp.max - start_date).days // self.ten_offsets_in_days)), 1000)

        def get_index_for_offset(off):
            self.start_date = Timestamp('1/1/1900')
            return date_range(self.start_date, periods=min(1000, get_period_count(self.start_date, off)), freq=off)
        self.idx = get_index_for_offset(MonthBegin(1, **{}))
        self.df = DataFrame(np.random.randn(len(self.idx), 10), index=self.idx)
        self.d = dict([(col, self.df[col]) for col in self.df.columns])

    def time_frame_ctor_dtindex_MonthBeginx1(self):
        DataFrame(self.d)


class frame_ctor_dtindex_MonthBeginx2(object):
    goal_time = 0.2

    def setup(self):

        def get_period_count(start_date, off):
            self.ten_offsets_in_days = ((start_date + (off * 10)) - start_date).days
            if (self.ten_offsets_in_days == 0):
                return 1000
            else:
                return min((9 * ((Timestamp.max - start_date).days // self.ten_offsets_in_days)), 1000)

        def get_index_for_offset(off):
            self.start_date = Timestamp('1/1/1900')
            return date_range(self.start_date, periods=min(1000, get_period_count(self.start_date, off)), freq=off)
        self.idx = get_index_for_offset(MonthBegin(2, **{}))
        self.df = DataFrame(np.random.randn(len(self.idx), 10), index=self.idx)
        self.d = dict([(col, self.df[col]) for col in self.df.columns])

    def time_frame_ctor_dtindex_MonthBeginx2(self):
        DataFrame(self.d)


class frame_ctor_dtindex_MonthEndx1(object):
    goal_time = 0.2

    def setup(self):

        def get_period_count(start_date, off):
            self.ten_offsets_in_days = ((start_date + (off * 10)) - start_date).days
            if (self.ten_offsets_in_days == 0):
                return 1000
            else:
                return min((9 * ((Timestamp.max - start_date).days // self.ten_offsets_in_days)), 1000)

        def get_index_for_offset(off):
            self.start_date = Timestamp('1/1/1900')
            return date_range(self.start_date, periods=min(1000, get_period_count(self.start_date, off)), freq=off)
        self.idx = get_index_for_offset(MonthEnd(1, **{}))
        self.df = DataFrame(np.random.randn(len(self.idx), 10), index=self.idx)
        self.d = dict([(col, self.df[col]) for col in self.df.columns])

    def time_frame_ctor_dtindex_MonthEndx1(self):
        DataFrame(self.d)


class frame_ctor_dtindex_MonthEndx2(object):
    goal_time = 0.2

    def setup(self):

        def get_period_count(start_date, off):
            self.ten_offsets_in_days = ((start_date + (off * 10)) - start_date).days
            if (self.ten_offsets_in_days == 0):
                return 1000
            else:
                return min((9 * ((Timestamp.max - start_date).days // self.ten_offsets_in_days)), 1000)

        def get_index_for_offset(off):
            self.start_date = Timestamp('1/1/1900')
            return date_range(self.start_date, periods=min(1000, get_period_count(self.start_date, off)), freq=off)
        self.idx = get_index_for_offset(MonthEnd(2, **{}))
        self.df = DataFrame(np.random.randn(len(self.idx), 10), index=self.idx)
        self.d = dict([(col, self.df[col]) for col in self.df.columns])

    def time_frame_ctor_dtindex_MonthEndx2(self):
        DataFrame(self.d)


class frame_ctor_dtindex_Nanox1(object):
    goal_time = 0.2

    def setup(self):

        def get_period_count(start_date, off):
            self.ten_offsets_in_days = ((start_date + (off * 10)) - start_date).days
            if (self.ten_offsets_in_days == 0):
                return 1000
            else:
                return min((9 * ((Timestamp.max - start_date).days // self.ten_offsets_in_days)), 1000)

        def get_index_for_offset(off):
            self.start_date = Timestamp('1/1/1900')
            return date_range(self.start_date, periods=min(1000, get_period_count(self.start_date, off)), freq=off)
        self.idx = get_index_for_offset(Nano(1, **{}))
        self.df = DataFrame(np.random.randn(len(self.idx), 10), index=self.idx)
        self.d = dict([(col, self.df[col]) for col in self.df.columns])

    def time_frame_ctor_dtindex_Nanox1(self):
        DataFrame(self.d)


class frame_ctor_dtindex_Nanox2(object):
    goal_time = 0.2

    def setup(self):

        def get_period_count(start_date, off):
            self.ten_offsets_in_days = ((start_date + (off * 10)) - start_date).days
            if (self.ten_offsets_in_days == 0):
                return 1000
            else:
                return min((9 * ((Timestamp.max - start_date).days // self.ten_offsets_in_days)), 1000)

        def get_index_for_offset(off):
            self.start_date = Timestamp('1/1/1900')
            return date_range(self.start_date, periods=min(1000, get_period_count(self.start_date, off)), freq=off)
        self.idx = get_index_for_offset(Nano(2, **{}))
        self.df = DataFrame(np.random.randn(len(self.idx), 10), index=self.idx)
        self.d = dict([(col, self.df[col]) for col in self.df.columns])

    def time_frame_ctor_dtindex_Nanox2(self):
        DataFrame(self.d)


class frame_ctor_dtindex_QuarterBeginx1(object):
    goal_time = 0.2

    def setup(self):

        def get_period_count(start_date, off):
            self.ten_offsets_in_days = ((start_date + (off * 10)) - start_date).days
            if (self.ten_offsets_in_days == 0):
                return 1000
            else:
                return min((9 * ((Timestamp.max - start_date).days // self.ten_offsets_in_days)), 1000)

        def get_index_for_offset(off):
            self.start_date = Timestamp('1/1/1900')
            return date_range(self.start_date, periods=min(1000, get_period_count(self.start_date, off)), freq=off)
        self.idx = get_index_for_offset(QuarterBegin(1, **{}))
        self.df = DataFrame(np.random.randn(len(self.idx), 10), index=self.idx)
        self.d = dict([(col, self.df[col]) for col in self.df.columns])

    def time_frame_ctor_dtindex_QuarterBeginx1(self):
        DataFrame(self.d)


class frame_ctor_dtindex_QuarterBeginx2(object):
    goal_time = 0.2

    def setup(self):

        def get_period_count(start_date, off):
            self.ten_offsets_in_days = ((start_date + (off * 10)) - start_date).days
            if (self.ten_offsets_in_days == 0):
                return 1000
            else:
                return min((9 * ((Timestamp.max - start_date).days // self.ten_offsets_in_days)), 1000)

        def get_index_for_offset(off):
            self.start_date = Timestamp('1/1/1900')
            return date_range(self.start_date, periods=min(1000, get_period_count(self.start_date, off)), freq=off)
        self.idx = get_index_for_offset(QuarterBegin(2, **{}))
        self.df = DataFrame(np.random.randn(len(self.idx), 10), index=self.idx)
        self.d = dict([(col, self.df[col]) for col in self.df.columns])

    def time_frame_ctor_dtindex_QuarterBeginx2(self):
        DataFrame(self.d)


class frame_ctor_dtindex_QuarterEndx1(object):
    goal_time = 0.2

    def setup(self):

        def get_period_count(start_date, off):
            self.ten_offsets_in_days = ((start_date + (off * 10)) - start_date).days
            if (self.ten_offsets_in_days == 0):
                return 1000
            else:
                return min((9 * ((Timestamp.max - start_date).days // self.ten_offsets_in_days)), 1000)

        def get_index_for_offset(off):
            self.start_date = Timestamp('1/1/1900')
            return date_range(self.start_date, periods=min(1000, get_period_count(self.start_date, off)), freq=off)
        self.idx = get_index_for_offset(QuarterEnd(1, **{}))
        self.df = DataFrame(np.random.randn(len(self.idx), 10), index=self.idx)
        self.d = dict([(col, self.df[col]) for col in self.df.columns])

    def time_frame_ctor_dtindex_QuarterEndx1(self):
        DataFrame(self.d)


class frame_ctor_dtindex_QuarterEndx2(object):
    goal_time = 0.2

    def setup(self):

        def get_period_count(start_date, off):
            self.ten_offsets_in_days = ((start_date + (off * 10)) - start_date).days
            if (self.ten_offsets_in_days == 0):
                return 1000
            else:
                return min((9 * ((Timestamp.max - start_date).days // self.ten_offsets_in_days)), 1000)

        def get_index_for_offset(off):
            self.start_date = Timestamp('1/1/1900')
            return date_range(self.start_date, periods=min(1000, get_period_count(self.start_date, off)), freq=off)
        self.idx = get_index_for_offset(QuarterEnd(2, **{}))
        self.df = DataFrame(np.random.randn(len(self.idx), 10), index=self.idx)
        self.d = dict([(col, self.df[col]) for col in self.df.columns])

    def time_frame_ctor_dtindex_QuarterEndx2(self):
        DataFrame(self.d)


class frame_ctor_dtindex_Secondx1(object):
    goal_time = 0.2

    def setup(self):

        def get_period_count(start_date, off):
            self.ten_offsets_in_days = ((start_date + (off * 10)) - start_date).days
            if (self.ten_offsets_in_days == 0):
                return 1000
            else:
                return min((9 * ((Timestamp.max - start_date).days // self.ten_offsets_in_days)), 1000)

        def get_index_for_offset(off):
            self.start_date = Timestamp('1/1/1900')
            return date_range(self.start_date, periods=min(1000, get_period_count(self.start_date, off)), freq=off)
        self.idx = get_index_for_offset(Second(1, **{}))
        self.df = DataFrame(np.random.randn(len(self.idx), 10), index=self.idx)
        self.d = dict([(col, self.df[col]) for col in self.df.columns])

    def time_frame_ctor_dtindex_Secondx1(self):
        DataFrame(self.d)


class frame_ctor_dtindex_Secondx2(object):
    goal_time = 0.2

    def setup(self):

        def get_period_count(start_date, off):
            self.ten_offsets_in_days = ((start_date + (off * 10)) - start_date).days
            if (self.ten_offsets_in_days == 0):
                return 1000
            else:
                return min((9 * ((Timestamp.max - start_date).days // self.ten_offsets_in_days)), 1000)

        def get_index_for_offset(off):
            self.start_date = Timestamp('1/1/1900')
            return date_range(self.start_date, periods=min(1000, get_period_count(self.start_date, off)), freq=off)
        self.idx = get_index_for_offset(Second(2, **{}))
        self.df = DataFrame(np.random.randn(len(self.idx), 10), index=self.idx)
        self.d = dict([(col, self.df[col]) for col in self.df.columns])

    def time_frame_ctor_dtindex_Secondx2(self):
        DataFrame(self.d)


class frame_ctor_dtindex_WeekOfMonthx1(object):
    goal_time = 0.2

    def setup(self):

        def get_period_count(start_date, off):
            self.ten_offsets_in_days = ((start_date + (off * 10)) - start_date).days
            if (self.ten_offsets_in_days == 0):
                return 1000
            else:
                return min((9 * ((Timestamp.max - start_date).days // self.ten_offsets_in_days)), 1000)

        def get_index_for_offset(off):
            self.start_date = Timestamp('1/1/1900')
            return date_range(self.start_date, periods=min(1000, get_period_count(self.start_date, off)), freq=off)
        self.idx = get_index_for_offset(WeekOfMonth(1, **{'week': 1, 'weekday': 1, }))
        self.df = DataFrame(np.random.randn(len(self.idx), 10), index=self.idx)
        self.d = dict([(col, self.df[col]) for col in self.df.columns])

    def time_frame_ctor_dtindex_WeekOfMonthx1(self):
        DataFrame(self.d)


class frame_ctor_dtindex_WeekOfMonthx2(object):
    goal_time = 0.2

    def setup(self):

        def get_period_count(start_date, off):
            self.ten_offsets_in_days = ((start_date + (off * 10)) - start_date).days
            if (self.ten_offsets_in_days == 0):
                return 1000
            else:
                return min((9 * ((Timestamp.max - start_date).days // self.ten_offsets_in_days)), 1000)

        def get_index_for_offset(off):
            self.start_date = Timestamp('1/1/1900')
            return date_range(self.start_date, periods=min(1000, get_period_count(self.start_date, off)), freq=off)
        self.idx = get_index_for_offset(WeekOfMonth(2, **{'week': 1, 'weekday': 1, }))
        self.df = DataFrame(np.random.randn(len(self.idx), 10), index=self.idx)
        self.d = dict([(col, self.df[col]) for col in self.df.columns])

    def time_frame_ctor_dtindex_WeekOfMonthx2(self):
        DataFrame(self.d)


class frame_ctor_dtindex_Weekx1(object):
    goal_time = 0.2

    def setup(self):

        def get_period_count(start_date, off):
            self.ten_offsets_in_days = ((start_date + (off * 10)) - start_date).days
            if (self.ten_offsets_in_days == 0):
                return 1000
            else:
                return min((9 * ((Timestamp.max - start_date).days // self.ten_offsets_in_days)), 1000)

        def get_index_for_offset(off):
            self.start_date = Timestamp('1/1/1900')
            return date_range(self.start_date, periods=min(1000, get_period_count(self.start_date, off)), freq=off)
        self.idx = get_index_for_offset(Week(1, **{}))
        self.df = DataFrame(np.random.randn(len(self.idx), 10), index=self.idx)
        self.d = dict([(col, self.df[col]) for col in self.df.columns])

    def time_frame_ctor_dtindex_Weekx1(self):
        DataFrame(self.d)


class frame_ctor_dtindex_Weekx2(object):
    goal_time = 0.2

    def setup(self):

        def get_period_count(start_date, off):
            self.ten_offsets_in_days = ((start_date + (off * 10)) - start_date).days
            if (self.ten_offsets_in_days == 0):
                return 1000
            else:
                return min((9 * ((Timestamp.max - start_date).days // self.ten_offsets_in_days)), 1000)

        def get_index_for_offset(off):
            self.start_date = Timestamp('1/1/1900')
            return date_range(self.start_date, periods=min(1000, get_period_count(self.start_date, off)), freq=off)
        self.idx = get_index_for_offset(Week(2, **{}))
        self.df = DataFrame(np.random.randn(len(self.idx), 10), index=self.idx)
        self.d = dict([(col, self.df[col]) for col in self.df.columns])

    def time_frame_ctor_dtindex_Weekx2(self):
        DataFrame(self.d)


class frame_ctor_dtindex_YearBeginx1(object):
    goal_time = 0.2

    def setup(self):

        def get_period_count(start_date, off):
            self.ten_offsets_in_days = ((start_date + (off * 10)) - start_date).days
            if (self.ten_offsets_in_days == 0):
                return 1000
            else:
                return min((9 * ((Timestamp.max - start_date).days // self.ten_offsets_in_days)), 1000)

        def get_index_for_offset(off):
            self.start_date = Timestamp('1/1/1900')
            return date_range(self.start_date, periods=min(1000, get_period_count(self.start_date, off)), freq=off)
        self.idx = get_index_for_offset(YearBegin(1, **{}))
        self.df = DataFrame(np.random.randn(len(self.idx), 10), index=self.idx)
        self.d = dict([(col, self.df[col]) for col in self.df.columns])

    def time_frame_ctor_dtindex_YearBeginx1(self):
        DataFrame(self.d)


class frame_ctor_dtindex_YearBeginx2(object):
    goal_time = 0.2

    def setup(self):

        def get_period_count(start_date, off):
            self.ten_offsets_in_days = ((start_date + (off * 10)) - start_date).days
            if (self.ten_offsets_in_days == 0):
                return 1000
            else:
                return min((9 * ((Timestamp.max - start_date).days // self.ten_offsets_in_days)), 1000)

        def get_index_for_offset(off):
            self.start_date = Timestamp('1/1/1900')
            return date_range(self.start_date, periods=min(1000, get_period_count(self.start_date, off)), freq=off)
        self.idx = get_index_for_offset(YearBegin(2, **{}))
        self.df = DataFrame(np.random.randn(len(self.idx), 10), index=self.idx)
        self.d = dict([(col, self.df[col]) for col in self.df.columns])

    def time_frame_ctor_dtindex_YearBeginx2(self):
        DataFrame(self.d)


class frame_ctor_dtindex_YearEndx1(object):
    goal_time = 0.2

    def setup(self):

        def get_period_count(start_date, off):
            self.ten_offsets_in_days = ((start_date + (off * 10)) - start_date).days
            if (self.ten_offsets_in_days == 0):
                return 1000
            else:
                return min((9 * ((Timestamp.max - start_date).days // self.ten_offsets_in_days)), 1000)

        def get_index_for_offset(off):
            self.start_date = Timestamp('1/1/1900')
            return date_range(self.start_date, periods=min(1000, get_period_count(self.start_date, off)), freq=off)
        self.idx = get_index_for_offset(YearEnd(1, **{}))
        self.df = DataFrame(np.random.randn(len(self.idx), 10), index=self.idx)
        self.d = dict([(col, self.df[col]) for col in self.df.columns])

    def time_frame_ctor_dtindex_YearEndx1(self):
        DataFrame(self.d)


class frame_ctor_dtindex_YearEndx2(object):
    goal_time = 0.2

    def setup(self):

        def get_period_count(start_date, off):
            self.ten_offsets_in_days = ((start_date + (off * 10)) - start_date).days
            if (self.ten_offsets_in_days == 0):
                return 1000
            else:
                return min((9 * ((Timestamp.max - start_date).days // self.ten_offsets_in_days)), 1000)

        def get_index_for_offset(off):
            self.start_date = Timestamp('1/1/1900')
            return date_range(self.start_date, periods=min(1000, get_period_count(self.start_date, off)), freq=off)
        self.idx = get_index_for_offset(YearEnd(2, **{}))
        self.df = DataFrame(np.random.randn(len(self.idx), 10), index=self.idx)
        self.d = dict([(col, self.df[col]) for col in self.df.columns])

    def time_frame_ctor_dtindex_YearEndx2(self):
        DataFrame(self.d)


class frame_ctor_list_of_dict(object):
    goal_time = 0.2

    def setup(self):
        (N, K) = (5000, 50)
        self.index = tm.makeStringIndex(N)
        self.columns = tm.makeStringIndex(K)
        self.frame = DataFrame(np.random.randn(N, K), index=self.index, columns=self.columns)
        try:
            self.data = self.frame.to_dict()
        except:
            self.data = self.frame.toDict()
        self.some_dict = self.data.values()[0]
        self.dict_list = [dict(zip(self.columns, row)) for row in self.frame.values]

    def time_frame_ctor_list_of_dict(self):
        DataFrame(self.dict_list)


class frame_ctor_nested_dict(object):
    goal_time = 0.2

    def setup(self):
        (N, K) = (5000, 50)
        self.index = tm.makeStringIndex(N)
        self.columns = tm.makeStringIndex(K)
        self.frame = DataFrame(np.random.randn(N, K), index=self.index, columns=self.columns)
        try:
            self.data = self.frame.to_dict()
        except:
            self.data = self.frame.toDict()
        self.some_dict = self.data.values()[0]
        self.dict_list = [dict(zip(self.columns, row)) for row in self.frame.values]

    def time_frame_ctor_nested_dict(self):
        DataFrame(self.data)


class frame_ctor_nested_dict_int64(object):
    goal_time = 0.2

    def setup(self):
        self.data = dict(((i, dict(((j, float(j)) for j in xrange(100)))) for i in xrange(2000)))

    def time_frame_ctor_nested_dict_int64(self):
        DataFrame(self.data)


class frame_from_series(object):
    goal_time = 0.2

    def setup(self):
        self.mi = MultiIndex.from_tuples([(x, y) for x in range(100) for y in range(100)])
        self.s = Series(randn(10000), index=self.mi)

    def time_frame_from_series(self):
        DataFrame(self.s)


class frame_get_numeric_data(object):
    goal_time = 0.2

    def setup(self):
        self.df = DataFrame(randn(10000, 25))
        self.df['foo'] = 'bar'
        self.df['bar'] = 'baz'
        self.df = self.df.consolidate()

    def time_frame_get_numeric_data(self):
        self.df._get_numeric_data()


class series_ctor_from_dict(object):
    goal_time = 0.2

    def setup(self):
        (N, K) = (5000, 50)
        self.index = tm.makeStringIndex(N)
        self.columns = tm.makeStringIndex(K)
        self.frame = DataFrame(np.random.randn(N, K), index=self.index, columns=self.columns)
        try:
            self.data = self.frame.to_dict()
        except:
            self.data = self.frame.toDict()
        self.some_dict = self.data.values()[0]
        self.dict_list = [dict(zip(self.columns, row)) for row in self.frame.values]

    def time_series_ctor_from_dict(self):
        Series(self.some_dict)