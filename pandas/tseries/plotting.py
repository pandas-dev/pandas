"""
Adapted from scikits.timeseries by Pierre GF Gerard-Marchant & Matt Knox
"""

#!!! TODO: Use the fact that axis can have units to simplify the process


from matplotlib import pylab
from matplotlib.axes import Subplot
from matplotlib.figure import Figure
from matplotlib.ticker import Formatter, Locator
from matplotlib.transforms import nonsingular

import numpy as np

import pandas.core.datetools as datetools

from pandas.core.datetools import Period
from pandas.core.index import PeriodIndex
from pandas.core.series import Series

import warnings

# Generic documentation ......................................................

_doc_parameters = dict(
figsize="""figsize : {None, tuple}
        Size of the figure, as a tuple (width, height) in inches.
        If None, defaults to rc figure.figsize.""",
dpi="""dpi : {None, int}, optional
        Resolution in dots per inches.
        If None, defaults to rc figure.dpi.""",
facecolor="""facecolor : {None, string}, optional
        Background color.
        If None, defaults to rc figure.facecolor.""",
edgecolor="""edgecolor : {None, string}, optional
        Border color.
        If None, defaults to rc figure.edgecolor.""",
linewidth="""linewidth : {float, None}
        Width of the patch edge line.""",
frameon="""frameon : {True, False}
        Whether to draw the frame around the figure.""",
subplotpars="""subplotpars : {None, var}
        A :class:`SubplotParams` instance, defaults to rc""",
mandatoryplotargs="""args : var
        Mandatory arguments for the creation of the subplot.
        These arguments should be given as ``nb_of_rows``, ``nb_of_columns``,
        ``plot_number``, or as a single 3-digit number if the 3 previous numbers
        are all lower than 10.""" )


#####---------------------------------------------------------------------------
#---- --- Matplotlib extensions ---
#####---------------------------------------------------------------------------

def add_generic_subplot(figure_instance, *args, **kwargs):
    """
    Generalizes the :meth:`matplotlib.Figure.add_subplot` method
    of :class:`~matplotlib.figure.Figure` to generic subplots.
    The specific Subplot object class to add is given through the keywords
    ``SubplotClass`` or ``class``.

    Parameters
    ----------
    figure_instance : Figure object
        Figure to which the generic subplot should be attached.
    args : {var}
        Miscellaneous arguments to the subplot.
    kwargs : {Dictionary}
        Optional keywords.
        The same keywords as ``Subplot`` are recognized, with the addition of:

        + *SubplotClass* : {string}
          Type of subplot.
        + *subclass* : {string}
          Shortcut to SubplotClass.
        + any keyword required by the ``SubplotClass`` subclass.

    """
    key = figure_instance._make_key(*args, ** kwargs)
    #!!!: Find why, sometimes, key is not hashable (even if tuple)
    # else, there's a fix below
    try:
        key.__hash__()
    except TypeError:
        key = str(key)

    if figure_instance._axstack.get(key):
        ax = figure_instance._axstack[key]
        figure_instance.sca(ax)
        return ax

    SubplotClass = kwargs.pop("SubplotClass", Subplot)
    SubplotClass = kwargs.pop("subclass", SubplotClass)
    if isinstance(args[0], Subplot):
        a = args[0]
        assert(a.get_figure() is figure_instance)
    else:
        a = SubplotClass(figure_instance, *args, **kwargs)

    figure_instance.axes.append(a)
    figure_instance._axstack.add(key, a)
    figure_instance.sca(a)
    return a

##### -------------------------------------------------------------------------
#---- --- Locators ---
##### -------------------------------------------------------------------------

def _get_default_annual_spacing(nyears):
    """
    Returns a default spacing between consecutive ticks for annual data.
    """
    if nyears < 11:
        (min_spacing, maj_spacing) = (1, 1)
    elif nyears < 20:
        (min_spacing, maj_spacing) = (1, 2)
    elif nyears < 50:
        (min_spacing, maj_spacing) = (1, 5)
    elif nyears < 100:
        (min_spacing, maj_spacing) = (5, 10)
    elif nyears < 200:
        (min_spacing, maj_spacing) = (5, 25)
    elif nyears < 600:
        (min_spacing, maj_spacing) = (10, 50)
    else:
        factor = nyears // 1000 + 1
        (min_spacing, maj_spacing) = (factor * 20, factor * 100)
    return (min_spacing, maj_spacing)


def period_break(dates, period):
    """
    Returns the indices where the given period changes.

    Parameters
    ----------
    dates : PeriodIndex
        Array of intervals to monitor.
    period : string
        Name of the period to monitor.
    """
    current = getattr(dates, period)
    previous = getattr(dates-1, period)
    return (current - previous).nonzero()[0]


def has_level_label(label_flags, vmin):
    """
    Returns true if the ``label_flags`` indicate there is at least one label
    for this level.

    if the minimum view limit is not an exact integer, then the first tick
    label won't be shown, so we must adjust for that.
    """
    if label_flags.size == 0 or (label_flags.size == 1 and
                                 label_flags[0] == 0 and
                                 vmin % 1 > 0.0):
        return False
    else:
        return True

FR_ANN = 1000
FR_QTR = 2000
FR_MTH = 3000
FR_WK = 4000
FR_BUS = 5000
FR_DAY = 6000
FR_HR = 7000
FR_MIN = 8000
FR_SEC = 9000
FR_UND = -10000

def get_freq_group(freq):
    if isinstance(freq, basestring):
        base, mult = datetools._get_freq_code(freq)
        freq = base
    return (freq // 1000) * 1000

def get_freq(freq):
    if isinstance(freq, basestring):
        base, mult = datetools._get_freq_code(freq)
        freq = base
    return freq

def _daily_finder(vmin, vmax, freq):
    periodsperday = -1

    if freq >= FR_HR:
        if freq == FR_SEC:
            periodsperday = 24 * 60 * 60
        elif freq == FR_MIN:
            periodsperday = 24 * 60
        elif freq == FR_HR:
            periodsperday = 24
        else:
            raise ValueError("unexpected frequency: %s" % freq)
        periodsperyear = 365 * periodsperday
        periodspermonth = 28 * periodsperday

    elif freq == FR_BUS:
        periodsperyear = 261
        periodspermonth = 19
    elif freq == FR_DAY:
        periodsperyear = 365
        periodspermonth = 28
    elif get_freq_group(freq) == FR_WK:
        periodsperyear = 52
        periodspermonth = 3
    elif freq == FR_UND:
        periodsperyear = 100
        periodspermonth = 10
    else:
        raise ValueError("unexpected frequency")

    # save this for later usage
    vmin_orig = vmin

    (vmin, vmax) = (Period(value=int(vmin), freq=freq),
                    Period(value=int(vmax), freq=freq))
    span = vmax.ordinal - vmin.ordinal + 1
    dates_ = PeriodIndex(start=vmin, end=vmax, freq=freq)
    # Initialize the output
    info = np.zeros(span,
                    dtype=[('val', int), ('maj', bool), ('min', bool),
                           ('fmt', '|S20')])
    info['val'][:] = dates_.values
    info['fmt'][:] = ''
    info['maj'][[0, -1]] = True
    # .. and set some shortcuts
    info_maj = info['maj']
    info_min = info['min']
    info_fmt = info['fmt']

    def first_label(label_flags):
        if (label_flags[0] == 0) and (label_flags.size > 1) and \
            ((vmin_orig % 1) > 0.0):
                return label_flags[1]
        else:
            return label_flags[0]

    # Case 1. Less than a month
    if span <= periodspermonth:

        day_start = period_break(dates_, 'day')
        month_start = period_break(dates_, 'month')

        def _hour_finder(label_interval, force_year_start):
            _hour = dates_.hour
            _prev_hour = (dates_-1).hour
            hour_start = (_hour - _prev_hour) != 0
            info_maj[day_start] = True
            info_min[hour_start & (_hour % label_interval == 0)] = True
            year_start = period_break(dates_, 'year')
            info_fmt[hour_start & (_hour % label_interval == 0)] = '%H:%M'
            info_fmt[day_start] = '%H:%M\n%d-%b'
            info_fmt[year_start] = '%H:%M\n%d-%b\n%Y'
            if force_year_start and not has_level_label(year_start, vmin_orig):
                info_fmt[first_label(day_start)] = '%H:%M\n%d-%b\n%Y'

        def _minute_finder(label_interval):
            hour_start = period_break(dates_, 'hour')
            _minute = dates_.minute
            _prev_minute = (dates_-1).minute
            minute_start = (_minute - _prev_minute) != 0
            info_maj[hour_start] = True
            info_min[minute_start & (_minute % label_interval == 0)] = True
            year_start = period_break(dates_, 'year')
            info_fmt = info['fmt']
            info_fmt[minute_start & (_minute % label_interval == 0)] = '%H:%M'
            info_fmt[day_start] = '%H:%M\n%d-%b'
            info_fmt[year_start] = '%H:%M\n%d-%b\n%Y'

        def _second_finder(label_interval):
            minute_start = period_break(dates_, 'minute')
            _second = dates_.second
            _prev_second = (dates_-1).second
            second_start = (_second - _prev_second) != 0
            info['maj'][minute_start] = True
            info['min'][second_start & (_second % label_interval == 0)] = True
            year_start = period_break(dates_, 'year')
            info_fmt = info['fmt']
            info_fmt[second_start & (_second % label_interval == 0)] = '%H:%M:%S'
            info_fmt[day_start] = '%H:%M:%S\n%d-%b'
            info_fmt[year_start] = '%H:%M:%S\n%d-%b\n%Y'

        if span < periodsperday / 12000.0: _second_finder(1)
        elif span < periodsperday / 6000.0: _second_finder(2)
        elif span < periodsperday / 2400.0: _second_finder(5)
        elif span < periodsperday / 1200.0: _second_finder(10)
        elif span < periodsperday / 800.0: _second_finder(15)
        elif span < periodsperday / 400.0: _second_finder(30)
        elif span < periodsperday / 150.0: _minute_finder(1)
        elif span < periodsperday / 70.0: _minute_finder(2)
        elif span < periodsperday / 24.0: _minute_finder(5)
        elif span < periodsperday / 12.0: _minute_finder(15)
        elif span < periodsperday / 6.0:  _minute_finder(30)
        elif span < periodsperday / 2.5: _hour_finder(1, False)
        elif span < periodsperday / 1.5: _hour_finder(2, False)
        elif span < periodsperday * 1.25: _hour_finder(3, False)
        elif span < periodsperday * 2.5: _hour_finder(6, True)
        elif span < periodsperday * 4: _hour_finder(12, True)
        else:
            info_maj[month_start] = True
            info_min[day_start] = True
            year_start = period_break(dates_, 'year')
            info_fmt = info['fmt']
            info_fmt[day_start] = '%d'
            info_fmt[month_start] = '%d\n%b'
            info_fmt[year_start] = '%d\n%b\n%Y'
            if not has_level_label(year_start, vmin_orig):
                if not has_level_label(month_start, vmin_orig):
                    info_fmt[first_label(day_start)] = '%d\n%b\n%Y'
                else:
                    info_fmt[first_label(month_start)] = '%d\n%b\n%Y'

    # Case 2. Less than three months
    elif span <= periodsperyear // 4:
        month_start = period_break(dates_, 'month')
        info_maj[month_start] = True
        if freq < FR_HR:
            info['min'] = True
        else:
            day_start = period_break(dates_, 'day')
            info['min'][day_start] = True
        week_start = period_break(dates_, 'week')
        year_start = period_break(dates_, 'year')
        info_fmt[week_start] = '%d'
        info_fmt[month_start] = '\n\n%b'
        info_fmt[year_start] = '\n\n%b\n%Y'
        if not has_level_label(year_start, vmin_orig):
            if not has_level_label(month_start, vmin_orig):
                info_fmt[first_label(week_start)] = '\n\n%b\n%Y'
            else:
                info_fmt[first_label(month_start)] = '\n\n%b\n%Y'
    # Case 3. Less than 14 months ...............
    elif span <= 1.15 * periodsperyear:
        year_start = period_break(dates_, 'year')
        month_start = period_break(dates_, 'month')
        week_start = period_break(dates_, 'week')
        info_maj[month_start] = True
        info_min[week_start] = True
        info_min[year_start] = False
        info_min[month_start] = False
        info_fmt[month_start] = '%b'
        info_fmt[year_start] = '%b\n%Y'
        if not has_level_label(year_start, vmin_orig):
            info_fmt[first_label(month_start)] = '%b\n%Y'
    # Case 4. Less than 2.5 years ...............
    elif span <= 2.5 * periodsperyear:
        year_start = period_break(dates_, 'year')
        quarter_start = period_break(dates_, 'quarter')
        month_start = period_break(dates_, 'month')
        info_maj[quarter_start] = True
        info_min[month_start] = True
        info_fmt[quarter_start] = '%b'
        info_fmt[year_start] = '%b\n%Y'
    # Case 4. Less than 4 years .................
    elif span <= 4 * periodsperyear:
        year_start = period_break(dates_, 'year')
        month_start = period_break(dates_, 'month')
        info_maj[year_start] = True
        info_min[month_start] = True
        info_min[year_start] = False

        month_break = dates_[month_start].month
        jan_or_jul = month_start[(month_break == 1) | (month_break == 7)]
        info_fmt[jan_or_jul] = '%b'
        info_fmt[year_start] = '%b\n%Y'
    # Case 5. Less than 11 years ................
    elif span <= 11 * periodsperyear:
        year_start = period_break(dates_, 'year')
        quarter_start = period_break(dates_, 'quarter')
        info_maj[year_start] = True
        info_min[quarter_start] = True
        info_min[year_start] = False
        info_fmt[year_start] = '%Y'
    # Case 6. More than 12 years ................
    else:
        year_start = period_break(dates_, 'year')
        year_break = dates_[year_start].years
        nyears = span / periodsperyear
        (min_anndef, maj_anndef) = _get_default_annual_spacing(nyears)
        major_idx = year_start[(year_break % maj_anndef == 0)]
        info_maj[major_idx] = True
        minor_idx = year_start[(year_break % min_anndef == 0)]
        info_min[minor_idx] = True
        info_fmt[major_idx] = '%Y'
    #............................................

    return info


def _monthly_finder(vmin, vmax, freq):
    if isinstance(freq, basestring):
        freq = get_freq(freq)

    if freq != FR_MTH:
        raise ValueError("Unexpected frequency")
    periodsperyear = 12

    vmin_orig = vmin
    (vmin, vmax) = (int(vmin), int(vmax))
    span = vmax - vmin + 1
    #..............
    # Initialize the output
    info = np.zeros(span,
                    dtype=[('val', int), ('maj', bool), ('min', bool),
                           ('fmt', '|S8')])
    info['val'] = np.arange(vmin, vmax + 1)
    dates_ = info['val']
    info['fmt'] = ''
    year_start = (dates_ % 12 == 1).nonzero()[0]
    info_maj = info['maj']
    info_fmt = info['fmt']
    #..............
    if span <= 1.15 * periodsperyear:
        info_maj[year_start] = True
        info['min'] = True

        info_fmt[:] = '%b'
        info_fmt[year_start] = '%b\n%Y'

        if not has_level_label(year_start, vmin_orig):
            if dates_.size > 1:
                idx = 1
            else:
                idx = 0
            info_fmt[idx] = '%b\n%Y'
    #..............
    elif span <= 2.5 * periodsperyear:
        quarter_start = (dates_ % 3 == 1).nonzero()
        info_maj[year_start] = True
        # TODO: Check the following : is it really info['fmt'] ?
        info['fmt'][quarter_start] = True
        info['min'] = True

        info_fmt[quarter_start] = '%b'
        info_fmt[year_start] = '%b\n%Y'
    #..............
    elif span <= 4 * periodsperyear:
        info_maj[year_start] = True
        info['min'] = True

        jan_or_jul = (dates_ % 12 == 1) | (dates_ % 12 == 7)
        info_fmt[jan_or_jul] = '%b'
        info_fmt[year_start] = '%b\n%Y'
    #..............
    elif span <= 11 * periodsperyear:
        quarter_start = (dates_ % 3 == 1).nonzero()
        info_maj[year_start] = True
        info['min'][quarter_start] = True

        info_fmt[year_start] = '%Y'
    #..................
    else:
        nyears = span / periodsperyear
        (min_anndef, maj_anndef) = _get_default_annual_spacing(nyears)
        years = dates_[year_start] // 12 + 1
        major_idx = year_start[(years % maj_anndef == 0)]
        info_maj[major_idx] = True
        info['min'][year_start[(years % min_anndef == 0)]] = True

        info_fmt[major_idx] = '%Y'
    #..............
    return info


def _quarterly_finder(vmin, vmax, freq):
    if isinstance(freq, basestring):
        freq = get_freq(freq)

    if get_freq_group(freq) != FR_QTR:
        raise ValueError("Unexpected frequency")
    periodsperyear = 4
    vmin_orig = vmin
    (vmin, vmax) = (int(vmin), int(vmax))
    span = vmax - vmin + 1
    #............................................
    info = np.zeros(span,
                    dtype=[('val', int), ('maj', bool), ('min', bool),
                           ('fmt', '|S8')])
    info['val'] = np.arange(vmin, vmax + 1)
    info['fmt'] = ''
    dates_ = info['val']
    info_maj = info['maj']
    info_fmt = info['fmt']
    year_start = (dates_ % 4 == 1).nonzero()[0]
    #..............
    if span <= 3.5 * periodsperyear:
        info_maj[year_start] = True
        info['min'] = True

        info_fmt[:] = 'Q%q'
        info_fmt[year_start] = 'Q%q\n%F'
        if not has_level_label(year_start, vmin_orig):
            if dates_.size > 1:
                idx = 1
            else:
                idx = 0
            info_fmt[idx] = 'Q%q\n%F'
    #..............
    elif span <= 11 * periodsperyear:
        info_maj[year_start] = True
        info['min'] = True
        info_fmt[year_start] = '%F'
    #..............
    else:
        years = dates_[year_start] // 4 + 1
        nyears = span / periodsperyear
        (min_anndef, maj_anndef) = _get_default_annual_spacing(nyears)
        major_idx = year_start[(years % maj_anndef == 0)]
        info_maj[major_idx] = True
        info['min'][year_start[(years % min_anndef == 0)]] = True
        info_fmt[major_idx] = '%F'
    #..............
    return info


def _annual_finder(vmin, vmax, freq):
    if isinstance(freq, basestring):
        freq = get_freq(freq)
    if get_freq_group(freq) != FR_ANN:
        raise ValueError("Unexpected frequency")
    (vmin, vmax) = (int(vmin), int(vmax + 1))
    span = vmax - vmin + 1
    #..............
    info = np.zeros(span,
                    dtype=[('val', int), ('maj', bool), ('min', bool),
                           ('fmt', '|S8')])
    info['val'] = np.arange(vmin, vmax + 1)
    info['fmt'] = ''
    dates_ = info['val']
    #..............
    (min_anndef, maj_anndef) = _get_default_annual_spacing(span)
    major_idx = dates_ % maj_anndef == 0
    info['maj'][major_idx] = True
    info['min'][(dates_ % min_anndef == 0)] = True
    info['fmt'][major_idx] = '%Y'
    #..............
    return info


def get_finder(freq):
    if isinstance(freq, basestring):
        freq = get_freq(freq)
    fgroup = get_freq_group(freq)

    if fgroup == FR_ANN:
        return _annual_finder
    elif fgroup == FR_QTR:
        return _quarterly_finder
    elif freq ==FR_MTH:
        return _monthly_finder
    elif (freq >= FR_BUS) or (freq == FR_UND) or fgroup == FR_WK:
        return _daily_finder
    else:
        errmsg = "Unsupported frequency: %s" % (freq)
        raise NotImplementedError(errmsg)

class TimeSeries_DateLocator(Locator):
    """
    Locates the ticks along an axis controlled by a :class:`Series`.

    Parameters
    ----------
    freq : {var}
        Valid frequency specifier.
    minor_locator : {False, True}, optional
        Whether the locator is for minor ticks (True) or not.
    dynamic_mode : {True, False}, optional
        Whether the locator should work in dynamic mode.
    base : {int}, optional
    quarter : {int}, optional
    month : {int}, optional
    day : {int}, optional
    """

    def __init__(self, freq, minor_locator=False, dynamic_mode=True,
                 base=1, quarter=1, month=1, day=1, plot_obj=None):
        if isinstance(freq, basestring):
            freq = get_freq(freq)
        self.freq = freq
        self.base = base
        (self.quarter, self.month, self.day) = (quarter, month, day)
        self.isminor = minor_locator
        self.isdynamic = dynamic_mode
        self.offset = 0
        self.plot_obj = plot_obj
        self.finder = get_finder(freq)

    def asminor(self):
        "Returns the locator set to minor mode."
        self.isminor = True
        return self

    def asmajor(self):
        "Returns the locator set to major mode."
        self.isminor = False
        return self

    def _get_default_locs(self, vmin, vmax):
        "Returns the default locations of ticks."

        if self.plot_obj.date_axis_info is None:
            self.plot_obj.date_axis_info = self.finder(vmin, vmax, self.freq)

        locator = self.plot_obj.date_axis_info

        if self.isminor:
            return np.compress(locator['min'], locator['val'])
        return np.compress(locator['maj'], locator['val'])

    def __call__(self):
        'Return the locations of the ticks.'

        vi = tuple(self.axis.get_view_interval())
        if vi != self.plot_obj.view_interval:
            self.plot_obj.date_axis_info = None
        self.plot_obj.view_interval = vi
        vmin, vmax = vi

        if vmax < vmin:
            vmin, vmax = vmax, vmin
        if self.isdynamic:
            locs = self._get_default_locs(vmin, vmax)
        else:
            base = self.base
            (d, m) = divmod(vmin, base)
            vmin = (d + 1) * base
            locs = range(vmin, vmax + 1, base)
        return locs

    def autoscale(self):
        """
        Sets the view limits to the nearest multiples of base that contain the
        data.
        """
        # requires matplotlib >= 0.98.0
        (vmin, vmax) = self.axis.get_data_interval()

        locs = self._get_default_locs(vmin, vmax)
        (vmin, vmax) = locs[[0, -1]]
        if vmin == vmax:
            vmin -= 1
            vmax += 1
        return nonsingular(vmin, vmax)

#####-------------------------------------------------------------------------
#---- --- Formatter ---
#####-------------------------------------------------------------------------
class TimeSeries_DateFormatter(Formatter):
    """
    Formats the ticks along an axis controlled by a :class:`PeriodIndex`.

    Parameters
    ----------
    freq : {int, string}
        Valid frequency specifier.
    minor_locator : {False, True}
        Whether the current formatter should apply to minor ticks (True) or
        major ticks (False).
    dynamic_mode : {True, False}
        Whether the formatter works in dynamic mode or not.
    """

    def __init__(self, freq, minor_locator=False, dynamic_mode=True,
                 plot_obj=None):
        if isinstance(freq, basestring):
            freq = get_freq(freq)
        self.format = None
        self.freq = freq
        self.locs = []
        self.formatdict = None
        self.isminor = minor_locator
        self.isdynamic = dynamic_mode
        self.offset = 0
        self.plot_obj = plot_obj
        self.finder = get_finder(freq)

    def asminor(self):
        "Returns the formatter set to minor mode."
        self.isminor = True
        return self

    def asmajor(self):
        "Returns the fromatter set to major mode."
        self.isminor = False
        return self

    def _set_default_format(self, vmin, vmax):
        "Returns the default ticks spacing."

        if self.plot_obj.date_axis_info is None:
            self.plot_obj.date_axis_info = self.finder(vmin, vmax, self.freq)
        info = self.plot_obj.date_axis_info

        if self.isminor:
            format = np.compress(info['min'] & np.logical_not(info['maj']),
                                 info)
        else:
            format = np.compress(info['maj'], info)
        self.formatdict = dict([(x, f) for (x, _, _, f) in format])
        return self.formatdict

    def set_locs(self, locs):
        'Sets the locations of the ticks'
        # don't actually use the locs. This is just needed to work with
        # matplotlib. Force to use vmin, vmax
        self.locs = locs

        (vmin, vmax) = vi = tuple(self.axis.get_view_interval())
        if vi != self.plot_obj.view_interval:
            self.plot_obj.date_axis_info = None
        self.plot_obj.view_interval = vi
        if vmax < vmin:
            (vmin, vmax) = (vmax, vmin)
        self._set_default_format(vmin, vmax)
    #
    def __call__(self, x, pos=0):
        if self.formatdict is None:
            return ''
        else:
            fmt = self.formatdict.pop(x, '')
            return Period(int(x), self.freq).strftime(fmt)


#####--------------------------------------------------------------------------
#---- --- TimeSeries plots ---
#####--------------------------------------------------------------------------

class TimeSeriesPlot(Subplot, object):
    """
    Based on : :class:`~matplotlib.axes.SubPlot`

    Defines a subclass of :class:`matplotlib.axes.Subplot` to plot time series.

    A :class:`~Series` is associated with the plot.  This time series is
    usually specified at the creation of the plot, through the optional
    parameter ``series``.  If no series is given at the creation, the first
    time series being plotted will be used as associated series.

    The associated series is stored in the :attr:`~TimeSeriesPlot.series`
    attribute.  It gives its frequency to the plot.  This frequency can be
    accessed through the attribute :attr:`freq`.  All the other series that
    will be plotted will be first converted to the :attr:`freq` frequency,
    using their :meth:`~asfreq` method.

    The same parameters used for the instanciation of a standard
    :class:`matplotlib.axes.Subplot` are recognized.

    Parameters
    ----------
    series : {None, Series}, optional
        The time series allocated to the plot.

    Attributes
    ----------
    freq : int
        Frequency of the plot.
    xdata : PeriodIndex
        The array of dates corresponding to the x axis.
    legendsymbols : list
    legendlabels : list
        List of the labels associated with each plot.
        The first label corresponds to the first plot, the second label to the
        second plot, and so forth.

    Warnings
    --------
    * Because the series to plot are first converted to the plot frequency, it
    * is recommended when plotting several series to associate the plot with
    * the series with the highest frequency, in order to keep a good level of
    * detail.

    """
    def __init__(self, fig=None, *args, **kwargs):

        # Retrieve the series ...................
        _series = kwargs.pop('series', getattr(fig, 'series', None))
        Subplot.__init__(self, fig, *args, **kwargs)

        # Process options .......................
        self.set_series(series=_series)

        self._austoscale = False
        # Get the data to plot
        self.legendsymbols = []
        self.legendlabels = []
        # keep track of axis format and tick info
        self.date_axis_info = None
        # used to keep track of current view interval to determine if we need
        # to reset date_axis_info
        self.view_interval = None


    def set_series(self, series=None):
        """
        Sets the time series associated with the plot.  If ``series`` is a
        valid :class:`~Series` object, the :attr:`xdata` attribute is updated
        to the ``_dates`` part of ``series``.
        """
        if series is not None:
            self._series = series.values
            self.xdata = series.index
        else:
            self._series = None
            self.xdata = None
    #
    def get_series(self):
        """
        Returns the data part of the time series associated with the plot,
        as a (subclass of) :class:`ndarray`.
        """
        return self._series
    #
    series = property(fget=get_series, fset=set_series,
                      doc="Underlying time series.")


    def set_ydata(self, series=None):
        errmsg = ("The use of 'set_ydata' is deprecated. "
                  "Please use 'set_series' instead")
        warnings.DepreciationWarning(errmsg)
        return self.set_series(series)
    #
    def get_ydata(self):
        errmsg = ("The use of 'get_ydata' is deprecated. "
                  "Please use 'get_series' instead")
        warnings.DepreciationWarning(errmsg)
        return self.get_series()
    #
    ydata = property(fget=get_ydata, fset=set_ydata,
                     doc="Underlying time series.")


    def get_freq(self):
        """
        Returns the underlying frequency of the plot
        """
        return getattr(self.xdata, 'freq', None)
    #
    freq = property(fget=get_freq, doc="Underlying frequency.")

    #......................................................
    def _check_plot_params(self, *args):
        """
        Defines the plot coordinates (and basic plotting arguments).
        """
        remaining = list(args)
        noinfo_msg = "No date information available!"
        # No args ? Use defaults, if any
        if len(args) == 0:
            if self.xdata is None:
                raise ValueError(noinfo_msg)
            return (self.xdata, self.series)
        output = []
        while len(remaining) > 0:
            a = remaining.pop(0)
            # The argument is a format: use default dates/
            if isinstance(a, str):
                if self.xdata is None:
                    raise ValueError(noinfo_msg)
                else:
                    output.extend([self.xdata, self.series, a])
            # The argument is a Series: use its dates for x
            elif isinstance(a, Series):
                (x, y) = (a.index, a.values)
                if len(remaining) > 0 and isinstance(remaining[0], str):
                    b = remaining.pop(0)
                    output.extend([x, y, b])
                else:
                    output.extend([x, y])
            # The argument is a PeriodIndex............
            elif isinstance(a, PeriodIndex):
                # Force to current freq
                if self.freq is not None:
                    if a.freq != self.freq:
                        a = a.asfreq(self.freq)
                # There's an argument after
                if len(remaining) > 0:
                    #...and it's a format string
                    if isinstance(remaining[0], str):
                        b = remaining.pop(0)
                        if self.series is None:
                            raise ValueError(noinfo_msg)
                        else:
                            output.extend([a, self.series, b])
                    #... and it's another date: use the default
                    elif isinstance(remaining[0], PeriodIndex):
                        if self.series is None:
                            raise ValueError(noinfo_msg)
                        else:
                            output.extend([a, self.series])
                    #... and it must be some data
                    else:
                        b = remaining.pop(0)
                        if len(remaining) > 0:
                            if isinstance(remaining[0], str):
                                c = remaining.pop(0)
                                output.extend([a, b, c])
                            else:
                                output.extend([a, b])
                else:
                    if self.series is None:
                        raise ValueError(noinfo_msg)
            # Otherwise..............................
            elif len(remaining) > 0 and isinstance(remaining[0], str):
                b = remaining.pop(0)
                if self.xdata is None:
                    raise ValueError(noinfo_msg)
                else:
                    output.extend([self.xdata, a, b])
            elif self.xdata is None:
                raise ValueError(noinfo_msg)
            else:
                output.extend([self.xdata, a])
        # Reinitialize the plot if needed ...........
        if self.xdata is None:
            self.xdata = output[0]
        # Force the xdata to the current frequency
        elif output[0].freq != self.freq:
            output = list(output)
            output[0] = output[0].asfreq(self.freq)
        return output
    #......................................................
    def tsplot(self, *args,  **kwargs):
        """
        Plots the data parsed in argument to the current axes.  This command
        accepts the same optional keywords as :func:`matplotlib.pyplot.plot`.

        The argument ``args`` is a variable length argument, allowing for
        multiple data to be plotted at once. Acceptable combinations are:

        No arguments or a format string: The time series associated with the
        subplot is plotted with the given format.  If no format string is
        given, the default format is used instead.  For example, to plot the
        underlying time series with the default format, use:

            >>> tsplot()

        To plot the underlying time series with a red solid line, use the
        command:

            >>> tsplot('r-')

        a :class:`~Series` object or one of its subclass with or without a
        format string: The given time series is plotted with the given format.
        If no format string is given, the default format is used instead.

        an array or sequence, with or without a format string: The data is
        plotted with the given format using the :attr:`~TimeSeriesPlot.xdata`
        attribute of the plot as abscissae.

        two arrays or sequences, with or without a format string: The data are
        plotted with the given format, using the first array as abscissae and
        the second as ordinates.


        Parameters
        ----------
        args : var
            Sequence of arguments, as described previously.
        kwargs : var
            Optional parameters.
            The same parameters are accepted as for
            :meth:`matplotlib.axes.Subplot.plot`.

        """
        args = self._check_plot_params(*args)
        self.legendlabels.append(kwargs.get('label', None))
        plotted = Subplot.plot(self, *args,  **kwargs)
        self.format_dateaxis()

        # when adding a right axis (using add_yaxis), for some reason the
        # x axis limits don't get properly set. This gets around the problem
        xlim = self.get_xlim()
        if xlim[0] == 0.0 and xlim[1] == 1.0:
            # if xlim still at default values, autoscale the axis
            self.autoscale_view()
        self.reset_datelimits()
        return plotted
    #......................................................
    def format_dateaxis(self):
        """
        Pretty-formats the date axis (x-axis).

        Major and minor ticks are automatically set for the frequency of the
        current underlying series.  As the dynamic mode is activated by
        default, changing the limits of the x axis will intelligently change
        the positions of the ticks.
        """
        # Get the locator class .................
        majlocator = TimeSeries_DateLocator(self.freq, dynamic_mode=True,
                                            minor_locator=False, plot_obj=self)
        minlocator = TimeSeries_DateLocator(self.freq, dynamic_mode=True,
                                            minor_locator=True, plot_obj=self)
        self.xaxis.set_major_locator(majlocator)
        self.xaxis.set_minor_locator(minlocator)
        # Get the formatter .....................
        majformatter = TimeSeries_DateFormatter(self.freq, dynamic_mode=True,
                                                minor_locator=False,
                                                plot_obj=self)
        minformatter = TimeSeries_DateFormatter(self.freq, dynamic_mode=True,
                                                minor_locator=True,
                                                plot_obj=self)
        self.xaxis.set_major_formatter(majformatter)
        self.xaxis.set_minor_formatter(minformatter)
        pylab.draw_if_interactive()
    #......................................................
    def set_dlim(self, start_date=None, end_date=None):
        """
        Sets the date limits of the plot to ``start_date`` and ``end_date``.
        The dates can be given as :class:`~Period` objects, strings or
        integers.

        Parameters
        ----------
        start_date : {var}
            Starting date of the plot. If None, the current left limit
            (earliest date) is used.
        end_date : {var}
            Ending date of the plot. If None, the current right limit (latest
            date) is used.
        """
        freq = self.freq
        if freq is None:
            raise ValueError("Undefined frequency! Date limits can't be set!")
        # TODO : Shouldn't we make get_datevalue a more generic function ?
        def get_datevalue(date, freq):
            if isinstance(date, Period):
                return date.asfreq(freq).value
            elif isinstance(date, str):
                return Period(date, freq).value
            elif isinstance(date, (int, float)) or \
                (isinstance(date, np.ndarray) and (date.size == 1)):
                return date
            elif date is None:
                return None
            raise ValueError("Unrecognizable date '%s'" % date)
        # Fix left limit ..............
        xleft = get_datevalue(start_date, freq)
        # Fix right limit .......
        xright = get_datevalue(end_date, freq)
        self.set_xlim(xleft, xright)
        return (xleft, xright)

    def reset_datelimits(self):
        """
        Reset the date range of the x axis to the date range of the underlying
        time series.
        """
        return self.set_xlim(self.xdata[[0, -1]])

    def get_dlim(self):
        """
        Returns the limits of the x axis as a :class:`~PeriodIndex`.
        """
        xlims = self.get_xlim()
        return PeriodIndex(xlims, freq=self.freq)

TSPlot = TimeSeriesPlot

def add_yaxis(fsp=None, position='right', yscale=None, basey=10, subsy=None):
    """
    Adds a second y-axis to a :class:`TimeSeriesPlot`.
    This function can also be used as a method.

    Parameters
    ----------
    fsp : {None, TimeSeriesPlot}
        Subplot to which the secondary y-axis is added.
        If None, the current subplot is selected: in that case, it should be
        a valid :class:`TimeSeriesPlot`.
        When used as a :class:`TimeSeriesPlot` method, this parameter points
        automatically to the calling subplot.
    position : {string}
        Position of the new axis, as either ``'left'`` or ``'right'``.
    yscale : {string}
        Scale of the new axis, as either ``'log'``, ``'linear'`` or ``None``.
        If None, uses the same scale as the first y axis.
    basey : {integer}
        Base of the logarithm for the new axis (if needed).
    subsy : {sequence}
        Sequence of the location of the minor ticks;
        None defaults to autosubs, which depend on the number of decades in
        the plot.
        Eg for base 10, ``subsy=(1,2,5)`` will  put minor ticks on 1, 2, 5, 11,
        12,15, 21, ....
        To turn off minor ticking, set ``subsy=[]``.

    Raises
    ------
    TypeError
        If the selected subplot is not a valid :class:`TimeSeriesPlot` object.

    """
    if fsp is None:
        fsp = pylab.gca()
    if not isinstance(fsp, TimeSeriesPlot):
        raise TypeError("The current plot is not a TimeSeriesPlot")
    fig = fsp.figure
    axisini = fsp.axis()
    fsp_alt_args = (fsp._rows, fsp._cols, fsp._num + 1)
    fsp_alt = fig.add_tsplot(frameon=False, position=fsp.get_position(),
                             sharex=fsp, *fsp_alt_args)
    # Set position ....................
    if position.lower() == 'right':
        (inipos, newpos) = ('left', 'right')
    else:
        (inipos, newpos) = ('right', 'left')
    # Force scales tics to one side ...
    fsp.yaxis.set_ticks_position(inipos)
    fsp.yaxis.set_label_position(inipos)
    # Force 2nd ticks to the other side..
    fsp_alt.yaxis.set_ticks_position(newpos)
    fsp_alt.yaxis.set_label_position(newpos)
    # Force period axis scale..........
    if yscale is None:
        yscale = fsp.get_yscale()
        try:
            basey = fsp.yaxis.get_major_locator()._base
        except AttributeError:
            basey = 10.
    fsp_alt.set_yscale(yscale, basey=basey, subsy=subsy)

    pylab.draw_if_interactive()
    return fsp_alt

TimeSeriesPlot.add_yaxis = add_yaxis

#####--------------------------------------------------------------------------
#---- --- TimeSeries Figures ---
#####--------------------------------------------------------------------------
class TimeSeriesFigure(Figure):
    """
    Based on :class:`matplotlib.figure.Figure`

    Create a new :class:`~matplotlib.figure.Figure` object.
    All the subplots share the same time series.

    The same parameters used for the creation of a standard
    :class:`~matplotlib.figure.Figure` are accepted.

    Parameters
    ----------
    series : {None, TimeSeries}, optional
        Underlying time series.
        All the subplots of the figure will share the same series.
    figsize : {None, tuple}
        Size of the figure, as a tuple (width, height) in inches.
        If None, defaults to rc figure.figsize.
    dpi : {None, int}, optional
        Resolution in dots per inches.
        If None, defaults to rc figure.dpi
    facecolor : {None, string}, optional
        Background color.
        If None, defaults to rc figure.facecolor.
    edgecolor : {None, string}, optional
        Border color.
        If None, defaults to rc figure.edgecolor.
    linewidth : {float, None}
        Width of the patch edge line.
    frameon : {True, False}
        Whether to draw the frame around the figure.

    """
    def __init__(self, **kwargs):
        self._series = series = kwargs.pop('series', None)
        Figure.__init__(self, **kwargs)
        fspnum = kwargs.pop('fspnum', None)
        if fspnum is not None:
            self.add_tsplot(fspnum, series=series)
    #.........
    def add_tsplot(self, *args, **kwargs):
        """
        Adds a :class:`TimeSeriesPlot` subplot to the current figure.

        Parameters
        ----------
        args : var
            Mandatory arguments for the creation of the subplot.
            These arguments should be given as ``nb_of_rows``,
            ``nb_of_columns``, ``plot_number``, or as a single 3-digit number
            if the 3 previous numbers are all lower than 10.
        kwargs : var
            Optional arguments, as recognized by `add_subplot`.
        """
        kwargs.update(SubplotClass=TimeSeriesPlot)
        if self._series is not None:
            kwargs.update(series=self._series)
        return add_generic_subplot(self, *args, **kwargs)

    add_subplot = add_tsplot

TSFigure = TimeSeriesFigure

#................................................
def tsfigure(num=None, figsize=None, dpi=None, facecolor=None,
             edgecolor=None, frameon=True, subplotpars=None,
             FigureClass=TSFigure):
    """
    Creates a new :class:`TimeSeriesFigure` object.

    Parameters
    ----------
    num : {None, int}, optional
        Number of the figure.
        If None, a new figure is created and ``num`` is incremented.
    %(figsize)s
    %(dpi)s
    %(facecolor)s
    %(edgecolor)s
    %(frameon)s
    %(subplotpars)s
    FigureClass : FigureClass
        Class of the figure to create
    """
    figargs = dict(num=num, figsize=figsize, dpi=dpi, facecolor=facecolor,
                   frameon=frameon, FigureClass=FigureClass,
                   subplotpars=subplotpars)
    fig = pylab.figure(**figargs)
    return fig

tsfigure.__doc__ %= _doc_parameters

def tsplot(series=None, num=None, figsize=None, dpi=None,
           facecolor=None, edgecolor=None, frameon=True, subplotpars=None,
           FigureClass=TSFigure):
    """
    Creates a new :class:`TimeSeriesFigure` object and plots a series

    Parameters
    ----------
    num : {None, int}, optional
        Number of the figure.
        If None, a new figure is created and ``num`` is incremented.
    %(figsize)s
    %(dpi)s
    %(facecolor)s
    %(edgecolor)s
    %(frameon)s
    %(subplotpars)s
    FigureClass : FigureClass
        Class of the figure to create

    Parameters
    ----------

    """
    # deal with kwargs
    fig = tsfigure(num=num, figsize=figsize, dpi=dpi,
                facecolor=facecolor, edgecolor=edgecolor,
                frameon=frameon, subplotpars=subplotpars,
                FigureClass=FigureClass)

    sub = fig.add_tsplot(111)
    ret = sub.tsplot(series)

    return ret

tsplot.__doc__ %= _doc_parameters

###############################################################################

