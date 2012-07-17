"""
Period formatters and locators adapted from scikits.timeseries by
Pierre GF Gerard-Marchant & Matt Knox
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
        freq, ax_freq, series = _maybe_resample(series, ax, freq, plotf,
                                                kwargs)

    # Convert DatetimeIndex to PeriodIndex
    if isinstance(series.index, DatetimeIndex):
        series = series.to_period(freq=freq)

    # Set ax with freq info
    _decorate_axes(ax, freq, kwargs)

    # mask missing values
    args = _maybe_mask(series)

    # how to make sure ax.clear() flows through?
    if not hasattr(ax, '_plot_data'):
        ax._plot_data = []
    ax._plot_data.append((series, kwargs))

    # styles
    style = kwargs.pop('style', None)
    if style is not None:
        args.append(style)

    lines = plotf(ax, *args,  **kwargs)
    label = kwargs.get('label', None)
    _reset_legend(ax, lines[0], label, kwargs)

    # set date formatter, locators and rescale limits
    format_dateaxis(ax, ax.freq)
    left, right = _get_xlim(ax.get_lines())
    ax.set_xlim(left, right)

    return lines

def _reset_legend(ax, line, label, kwargs):
    ax, leg = _get_ax_legend(ax)
    if leg and (kwargs.get('legend', True)):
        ext_lines = leg.get_lines()
        ext_labels = [x.get_text() for x in leg.get_texts()]
        title = leg.get_title().get_text()
        if title == 'None':
            title = None

        ext_lines.append(line)
        ext_labels.append(label)
        ax.legend(ext_lines, ext_labels, loc='best', title=title)

def _get_ax_legend(ax):
    leg = ax.get_legend()

    other_ax = getattr(ax, 'right_ax', None) or getattr(ax, 'left_ax', None)
    other_leg = None
    if other_ax is not None:
        other_leg = other_ax.get_legend()

    if leg is None:
        leg = other_leg
        ax = other_ax

    return ax, leg

def _maybe_resample(series, ax, freq, plotf, kwargs):
    ax_freq = getattr(ax, 'freq', None)
    if (ax_freq is not None) and (freq != ax_freq):
        if frequencies.is_subperiod(freq, ax_freq): # upsample existing
            _upsample_others(ax, freq, ax_freq, plotf, kwargs)
            ax_freq = freq
        elif frequencies.is_superperiod(freq, ax_freq): # upsample input
            series = series.asfreq(ax_freq).dropna()
            freq = ax_freq
        elif _is_sup(freq, ax_freq): # one is weekly
            how = kwargs.pop('how', 'last')
            series = series.resample('D', how=how).dropna()
            series = series.resample(ax_freq, how=how).dropna()
            freq = ax_freq
        elif _is_sub(freq, ax_freq):
            _upsample_others(ax, freq, ax_freq, plotf, kwargs, True)
            ax_freq = freq
        else:
            raise ValueError('Incompatible frequency conversion')
    return freq, ax_freq, series

def _is_sub(f1, f2):
    return ((f1.startswith('W') and frequencies.is_subperiod('D', f2)) or
            (f2.startswith('W') and frequencies.is_subperiod(f1, 'D')))

def _is_sup(f1, f2):
    return ((f1.startswith('W') and frequencies.is_superperiod('D', f2)) or
            (f2.startswith('W') and frequencies.is_superperiod(f1, 'D')))

def _upsample_others(ax, freq, ax_freq, plotf, kwargs,
                     via_daily=False):
    legend = ax.get_legend()
    lines, labels = _replot_ax(ax, freq, ax_freq, plotf, kwargs, via_daily)

    other_ax = None
    if hasattr(ax, 'left_ax'):
        other_ax = ax.left_ax
    if hasattr(ax, 'right_ax'):
        other_ax = ax.right_ax

    if other_ax is not None:
        other_leg = other_ax.get_legend()
        rlines, rlabels = _replot_ax(other_ax, freq, ax_freq, plotf, kwargs,
                                     via_daily)
        lines.extend(rlines)
        labels.extend(rlabels)

    if (legend is not None and kwargs.get('legend', True) and
        len(lines) > 0):
        title = legend.get_title().get_text()
        if title == 'None':
            title = None
        ax.legend(lines, labels, loc='best', title=title)

def _replot_ax(ax, freq, ax_freq, plotf, kwargs, via_daily):
    data = ax._plot_data
    ax._plot_data = []
    ax.clear()
    _decorate_axes(ax, freq, kwargs)
    lines = []
    labels = []
    for series, kwds in data:
        series = _upsample(series, freq, via_daily)
        ax._plot_data.append(series)
        args = _maybe_mask(series)
        lines.append(plotf(ax, *args, **kwds)[0])
        labels.append(com._stringify(series.name))

    return lines, labels

def _upsample(series, freq, via_daily):
    if not via_daily:
        return series.resample(freq).dropna()
    else:
        return series.resample('D').resample(freq).dropna()

def _decorate_axes(ax, freq, kwargs):
    ax.freq = freq
    xaxis = ax.get_xaxis()
    xaxis.freq = freq
    if not hasattr(ax, 'legendlabels'):
        ax.legendlabels = [kwargs.get('label', None)]
    else:
        ax.legendlabels.append(kwargs.get('label', None))
    ax.view_interval = None
    ax.date_axis_info = None

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
