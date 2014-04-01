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
        info = ''
        if self.year is not None:
            info += 'year=%s, ' % self.year
        info += 'month=%s, day=%s, ' % (self.month, self.day)
        
        if self.offset is not None:
            info += 'offset=%s' % self.offset
            
        if self.observance is not None:
            info += 'observance=%s' % self.observance
        
        repr = 'Holiday: %s (%s)' % (self.name, info)
        return repr
    
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

holiday_calendars = {}
def register(cls):
    try:
        name = cls.name
    except:
        name = cls.__name__
    holiday_calendars[name] = cls
    
def get_calendar(name):
    '''
    Return an instance of a calendar based on its name.
    
    Parameters
    ----------
    name : str
        Calendar name to return an instance of
    '''
    return holiday_calendars[name]()

class HolidayMetaClass(type):
    def __new__(cls, clsname, bases, attrs):
        calendar_class = super(HolidayMetaClass, cls).__new__(cls, clsname, bases, attrs)
        register(calendar_class)
        return calendar_class

class AbstractHolidayCalendar(object):
    '''
    Abstract interface to create holidays following certain rules.
    '''
    __metaclass__ = HolidayMetaClass
    rules = []
    
    def __init__(self, name=None, rules=None):
        '''
        Initializes holiday object with a given set a rules.  Normally
        classes just have the rules defined within them.
        
        Parameters
        ----------
        name : str 
            Name of the holiday calendar, defaults to class name
        rules : array of Holiday objects
            A set of rules used to create the holidays.
        '''
        super(AbstractHolidayCalendar, self).__init__()
        if name is None:
            name = self.__class__.__name__
        self.name = name
        
        if rules is not None:
            self.rules = rules

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
            
        if self.rules is None:
            raise Exception('Holiday Calendar %s does not have any '\
                            'rules specified' % self.name)
        
        if return_names:
            holidays = None
        else:
            holidays    = []
        for rule in self.rules:
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

    @staticmethod
    def merge_class(base, other):
        '''
        Merge holiday calendars together.  The base calendar
        will take precedence to other. The merge will be done
        based on each holiday's name.
        
        Parameters
        ----------
        base : AbstractHolidayCalendar instance of array of Holiday objects
        other : AbstractHolidayCalendar instance or array of Holiday objects
        '''
        if isinstance(other, AbstractHolidayCalendar):
            other = other.rules
        if not isinstance(other, list):
            other = [other]
        other_holidays = {holiday.name: holiday for holiday in other}
            
        if isinstance(base, AbstractHolidayCalendar):
            base = base.rules
        if not isinstance(base, list):
            base = [base]
        base_holidays = {holiday.name: holiday for holiday in base}
        
        other_holidays.update(base_holidays)
        return other_holidays.values()

    def merge(self, other, inplace=False):
        '''
        Merge holiday calendars together.  The caller's class
        rules take precedence.  The merge will be done
        based on each holiday's name.
        
        Parameters
        ----------
        other : holiday calendar
        inplace : bool (default=False)
            If True set rule_table to holidays, else return array of Holidays
        '''
        holidays    =   self.merge_class(self, other)
        if inplace:
            self.rules = holidays
        else:
            return holidays

USMemorialDay     = Holiday('MemorialDay', month=5, day=24, 
                            offset=DateOffset(weekday=MO(1)))
USLaborDay        = Holiday('Labor Day', month=9, day=1, 
                            offset=DateOffset(weekday=MO(1)))
USColumbusDay     = Holiday('Columbus Day', month=10, day=1,
                            offset=DateOffset(weekday=MO(2)))
USThanksgivingDay = Holiday('Thanksgiving', month=11, day=1, 
                            offset=DateOffset(weekday=TH(4)))
USMartinLutherKingJr = Holiday('Dr. Martin Luther King Jr.', month=1, day=1, 
                               offset=DateOffset(weekday=MO(3)))
USPresidentsDay      = Holiday('President''s Day', month=2, day=1, 
                               offset=DateOffset(weekday=MO(3)))

class USFederalHolidayCalendar(AbstractHolidayCalendar):
    '''
    US Federal Government Holiday Calendar based on rules specified
    by: https://www.opm.gov/policy-data-oversight/snow-dismissal-procedures/federal-holidays/
    '''
    rules = [ 
        Holiday('New Years Day', month=1,  day=1,  observance=Nearest), 
        USMartinLutherKingJr,
        USPresidentsDay,
        USMemorialDay,
        Holiday('July 4th', month=7,  day=4,  observance=Nearest),
        USLaborDay,
        USColumbusDay,
        Holiday('Veterans Day', month=11, day=11, observance=Nearest),
        USThanksgivingDay,
        Holiday('Christmas', month=12, day=25, observance=Nearest)
        ]

class NERCHolidayCalendar(AbstractHolidayCalendar):
    '''
    North American Electric Reliability Corporation holiday rules:
    http://www.nerc.com/comm/OC/RS%20Agendas%20Highlights%20and%20Minutes%20DL/Additional_Off-peak_Days.pdf
    '''
    rules = [
        Holiday('New Years Day', month=1, day=1, observance=Sunday),
        USMemorialDay,
        Holiday('July 4th', month=7,  day=4, observance=Sunday),
        USLaborDay,
        USThanksgivingDay,
        Holiday('Christmas', month=12, day=25, observance=Sunday)
        ]
    
def HolidayCalendarFactory(name, base, other, base_class=AbstractHolidayCalendar):
    rules = AbstractHolidayCalendar.merge_class(base, other)
    calendar_class = type(name, (base_class,), {"rules": rules, "name": name})
    return calendar_class
