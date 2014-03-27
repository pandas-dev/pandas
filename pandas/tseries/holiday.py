from pandas import DateOffset, date_range, DatetimeIndex, Series
from datetime import datetime
from pandas.tseries.offsets import Easter
from dateutil.relativedelta import MO, TU, WE, TH, FR, SA, SU

def Sunday(dt):
    '''
    If the holiday falls on Sunday, make Monday a holiday (nothing
    happens for Saturday.
    '''
    if dt.isoweekday() == 7:
        return dt + DateOffset(+1)
    else:
        return dt
    
def Nearest(dt):
    '''
    If the holiday falls on a weekend, make it a 3-day weekend by making
    Saturday a Friday holiday and Sunday a Monday holiday.
    '''
    if dt.isoweekday() == 6:
        return dt + DateOffset(-1)
    elif dt.isoweekday() == 7:
        return dt + DateOffset(+1)
    else:
        return dt

#TODO: Need to add an observance function when a holiday
# falls on a Tuesday and get a 4-day weekend
# def Nearest4(dt):
#     '''
#     If the holiday falls on Tuesday,
#     make Monday a holiday as well, otherwise
#     follow the rules for Nearest (a
#     3-day weekend).
#     '''
#     if dt.isoweekday() == 2:
#         return dt - DateOffset()
#     else:
#         return Nearest(dt)

class Holiday(object):
    '''
    Class that defines a holiday with start/end dates and rules
    for observance.
    '''
    def __init__(self, name, year=None, month=None, day=None, offset=None,
                 observance=None, start_date=None, end_date=None):
        self.name   =   name
        self.year   =   year
        self.month  =   month
        self.day    =   day
        self.offset =   offset
        self.start_date = start_date
        self.end_date   = end_date
        self.observance = observance
    
    def __repr__(self):
        #FIXME: This should handle observance rules as well
        return 'Holiday %s (%s, %s, %s)' % (self.name, self.month, self.day, 
                                            self.offset)
    
    def dates(self, start_date, end_date):
        
        if self.year is not None:
            return datetime(self.year, self.month, self.day)
        
        if self.start_date is not None:
            start_date = self.start_date
            
        if self.end_date is not None:
            end_date = self.end_date
        
        year_offset = DateOffset(years=1)   
        baseDate = datetime(start_date.year, self.month, self.day)
        dates = date_range(baseDate, end_date, freq=year_offset)
        
        return self._apply_rule(dates)
    
    def dates_with_name(self, start_date, end_date):
        
        dates = self.dates(start_date, end_date)
        return Series(self.name, index=dates)

    def _apply_rule(self, dates):   
        '''
        Apply the given offset/observance to an 
        iterable of dates.
        
        Parameters
        ----------
        dates : array-like
            Dates to apply the given offset/observance rule
        
        Returns
        -------
        Dates with rules applied
        '''
        if self.observance is not None:
            return map(lambda d: self.observance(d), dates)

        if not isinstance(self.offset, list):
            offsets =   [self.offset]
        else:
            offsets =   self.offset
            
        for offset in offsets:
            dates = map(lambda d: d + offset, dates)
            
        return dates

class AbstractHolidayCalendar(object):
    '''
    Abstract interface to create holidays following certain rules.
    '''
    _rule_table = []
    
    def __init__(self, rules=None):
        '''
        Initializes holiday object with a given set a rules.  Normally
        classes just have the rules defined within them.
        
        Parameters
        ----------
        rules : array of Holiday objects
            A set of rules used to create the holidays.
        '''
        super(AbstractHolidayCalendar, self).__init__()
        if rules is not None:
            self._rule_table = rules
        
    @property
    def holiday_rules(self):
        return self._rule_table

    def holidays(self, start=None, end=None, return_names=False):
        '''
        Returns a curve with holidays between start_date and end_date
        
        Parameters
        ----------
        start : starting date, datetime-like, optional
        end : ending date, datetime-like, optional
        return_names : bool, optional
            If True, return a series that has dates and holiday names.
            False will only return a DatetimeIndex of dates.

        Returns
        -------
            DatetimeIndex of holidays
        '''
        #FIXME: Where should the default limits exist?
        if start is None:
            start = datetime(1970, 1, 1)
            
        if end is None:
            end = datetime(2030, 12, 31)
            
        if self.holiday_rules is None:
            raise Exception('Holiday Calendar %s does not have any '\
                            'rules specified' % self.calendarName)
        
        if return_names:
            holidays = None
        else:
            holidays    = []
        for rule in self.holiday_rules:
            if return_names:
                rule_holidays = rule.dates_with_name(start, end)
                if holidays is None:
                    holidays = rule_holidays
                else:
                    holidays = holidays.append(rule_holidays)
            else:
                holidays += rule.dates(start, end)

        if return_names:
            return holidays.sort_index()
        else:
            return DatetimeIndex(holidays).order(False)

USMemorialDay     = Holiday('MemorialDay', month=5, day=24, 
                            offset=DateOffset(weekday=MO(1)))
USLaborDay        = Holiday('Labor Day', month=9, day=1, 
                            offset=DateOffset(weekday=MO(1)))
USThanksgivingDay = Holiday('Thanksgiving', month=11, day=1, 
                            offset=DateOffset(weekday=TH(4)))
USMartinLutherKingJr = Holiday('Dr. Martin Luther King Jr.', month=1, day=1, 
                               offset=DateOffset(weekday=MO(3)))
USPresidentsDay      = Holiday('President''s Day', month=2, day=1, 
                               offset=DateOffset(weekday=MO(3)))

class USFederalHolidayCalendar(AbstractHolidayCalendar):
    
    _rule_table = [ 
        Holiday('New Years Day', month=1,  day=1,  observance=Nearest), 
        USMartinLutherKingJr,
        USPresidentsDay,
        USMemorialDay,
        Holiday('July 4th', month=7,  day=4,  observance=Nearest),
        USLaborDay, 
        Holiday('Columbus Day', month=10, day=1,  offset=DateOffset(weekday=MO(2))),
        Holiday('Veterans Day', month=11, day=11, observance=Nearest),
        USThanksgivingDay,
        Holiday('Christmas', month=12, day=25, observance=Nearest)
        ]

class NERCHolidayCalendar(AbstractHolidayCalendar):
    
    _rule_table = [
        Holiday('New Years Day', month=1, day=1, observance=Sunday),
        USMemorialDay,
        Holiday('July 4th', month=7,  day=4, observance=Sunday),
        USLaborDay,
        USThanksgivingDay,
        Holiday('Christmas', month=12, day=25, observance=Sunday)
        ]
