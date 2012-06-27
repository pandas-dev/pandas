"""
Adapted from scikits.timeseries by Pierre GF Gerard-Marchant & Matt Knox
"""

#!!! TODO: Use the fact that axis can have units to simplify the process
import datetime as pydt
from datetime import datetime

from matplotlib import pylab
import matplotlib.units as units

import numpy as np

from pandas import isnull
from pandas.tseries.period import Period
from pandas.tseries.offsets import DateOffset
import pandas.tseries.frequencies as frequencies
from pandas.tseries.index import DatetimeIndex
import pandas.core.common as com

from pandas.tseries.converter import (PeriodConverter, TimeSeries_DateLocator,
                                      TimeSeries_DateFormatter)

units.registry[Period] = PeriodConverter()
#----------------------------------------------------------------------
# Plotting functions and monkey patches

def tsplot(series, plotf, **kwargs):
    """
    Plots a Series on the given Matplotlib axes or the current axes

    Parameters
    ----------
    axes : Axes
    series : Series

    Notes
    _____
    Supports same kwargs as Axes.plot

    """
    # Used inferred freq is possible, need a test case for inferred
    if 'ax' in kwargs:
        ax = kwargs.pop('ax')
    else:
        import matplotlib.pyplot as plt
        ax = plt.gca()

    freq = _get_freq(ax, series)
    # resample against axes freq if necessary
    if freq is None: # pragma: no cover
        raise ValueError('Cannot use dynamic axis without frequency info')
    else:
        ax_freq = getattr(ax, 'freq', None)
        if (ax_freq is not None) and (freq != ax_freq):
            if frequencies.is_subperiod(freq, ax_freq): # downsample
                how = kwargs.pop('how', 'last')
                series = series.resample(ax_freq, how=how)
            elif frequencies.is_superperiod(freq, ax_freq):
                series = series.resample(ax_freq)
            else: # one freq is weekly
                how = kwargs.pop('how', 'last')
                series = series.resample('D', how=how, fill_method='pad')
                series = series.resample(ax_freq, how=how, fill_method='pad')
            freq = ax_freq

    # Convert DatetimeIndex to PeriodIndex
    if isinstance(series.index, DatetimeIndex):
        series = series.to_period(freq=freq)

    style = kwargs.pop('style', None)

    # Specialized ts plotting attributes for Axes
    ax.freq = freq
    xaxis = ax.get_xaxis()
    xaxis.freq = freq
    ax.legendlabels = [kwargs.get('label', None)]
    ax.view_interval = None
    ax.date_axis_info = None

    # format args and lot
    args = _maybe_mask(series)

    if style is not None:
        args.append(style)

    plotf(ax, *args,  **kwargs)

    format_dateaxis(ax, ax.freq)

    left, right = _get_xlim(ax.get_lines())
    ax.set_xlim(left, right)

    return ax

def _maybe_mask(series):
    mask = isnull(series)
    if mask.any():
        masked_array = np.ma.array(series.values)
        masked_array = np.ma.masked_where(mask, masked_array)
        args = [series.index, masked_array]
    else:
        args = [series.index, series]
    return args

def _get_freq(ax, series):
    # get frequency from data
    freq = getattr(series.index, 'freq', None)
    if freq is None:
        freq = getattr(series.index, 'inferred_freq', None)

    ax_freq = getattr(ax, 'freq', None)

    # use axes freq if no data freq
    if freq is None:
        freq = ax_freq

    # get the period frequency
    if isinstance(freq, DateOffset):
        freq = freq.rule_code
    else:
        freq = frequencies.get_base_alias(freq)

    freq = frequencies.get_period_alias(freq)

    return freq

def _get_xlim(lines):
    left, right = np.inf, -np.inf
    for l in lines:
        x = l.get_xdata()
        left = min(x[0].ordinal, left)
        right = max(x[-1].ordinal, right)
    return left, right

def get_datevalue(date, freq):
    if isinstance(date, Period):
        return date.asfreq(freq).ordinal
    elif isinstance(date, (str, datetime, pydt.date, pydt.time)):
        return Period(date, freq).ordinal
    elif (com.is_integer(date) or com.is_float(date) or
          (isinstance(date, np.ndarray) and (date.size == 1))):
        return date
    elif date is None:
        return None
    raise ValueError("Unrecognizable date '%s'" % date)

# Patch methods for subplot. Only format_dateaxis is currently used.
# Do we need the rest for convenience?

def format_dateaxis(subplot, freq):
    """
    Pretty-formats the date axis (x-axis).

    Major and minor ticks are automatically set for the frequency of the
    current underlying series.  As the dynamic mode is activated by
    default, changing the limits of the x axis will intelligently change
    the positions of the ticks.
    """
    majlocator = TimeSeries_DateLocator(freq, dynamic_mode=True,
                                        minor_locator=False,
                                        plot_obj=subplot)
    minlocator = TimeSeries_DateLocator(freq, dynamic_mode=True,
                                        minor_locator=True,
                                        plot_obj=subplot)
    subplot.xaxis.set_major_locator(majlocator)
    subplot.xaxis.set_minor_locator(minlocator)

    majformatter = TimeSeries_DateFormatter(freq, dynamic_mode=True,
                                            minor_locator=False,
                                            plot_obj=subplot)
    minformatter = TimeSeries_DateFormatter(freq, dynamic_mode=True,
                                            minor_locator=True,
                                            plot_obj=subplot)
    subplot.xaxis.set_major_formatter(majformatter)
    subplot.xaxis.set_minor_formatter(minformatter)
    pylab.draw_if_interactive()




