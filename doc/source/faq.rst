.. currentmodule:: pandas
.. _faq:

********************************
Frequently Asked Questions (FAQ)
********************************

Migrating from scikits.timeseries to pandas >= 0.8.0
----------------------------------------------------

.. csv-table::
    :header: "scikits.timeseries", "pandas", "Notes"
    :widths: 20, 20, 60

    Date, Period,
    DateArray, PeriodIndex,
    convert, resample,
    convert_to_annual, pivot_annual,
