.. currentmodule:: pandas
.. _faq:

********************************
Frequently Asked Questions (FAQ)
********************************

Migrating from scikits.timeseries to pandas >= 0.8.0
----------------------------------------------------

Starting with pandas 0.8.0, users of scikits.timeseries should have all of the
features that they need to migrate their code to use pandas. Portions of the
scikits.timeseries codebase for implementing calendar logic and timespan
frequency conversions (but **not** resampling, that has all been implemented
from scratch from the ground up) have been ported to the pandas codebase.

The scikits.timeseries notions of ``Date`` and ``DateArray`` are responsible
for implementing calendar logic:

::

    In [16]: dt = ts.Date('Q', '1984Q3')

    # sic
    In [17]: dt
    Out[17]: <Q-DEC : 1984Q1>

    In [18]: dt.asfreq('D', 'start')
    Out[18]: <D : 01-Jan-1984>

    In [19]: dt.asfreq('D', 'end')
    Out[19]: <D : 31-Mar-1984>

    In [20]: dt + 3
    Out[20]: <Q-DEC : 1984Q4>

``Date`` and ``DateArray`` from scikits.timeseries have been reincarnated in
pandas ``Period`` and ``PeriodIndex``:

.. ipython:: python

    p = Period('1984Q3')
    p
    p.asfreq('D', 'start')
    p.asfreq('D', 'end')
    (p + 3).asfreq('T') + 6 * 60 + 30
    rng = period_range('1990', '2010', freq='A')
    rng
    rng.asfreq('B', 'end') - 3

.. csv-table::
    :header: "scikits.timeseries", "pandas", "Notes"
    :widths: 20, 20, 60

    Date, Period, "A span of time, from yearly through to secondly"
    DateArray, PeriodIndex, "An array of timespans"
    convert, resample, "Frequency conversion in scikits.timeseries"
    convert_to_annual, pivot_annual,

Converting to and from period format
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Resampling with timestamps and periods
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
