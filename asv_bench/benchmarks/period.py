from pandas import PeriodIndex, date_range


class create_period_index_from_date_range(object):
    goal_time = 0.2

    def time_period_index(self):
        # Simulate irregular PeriodIndex
        PeriodIndex(date_range('1985', periods=1000).to_pydatetime(), freq='D')
