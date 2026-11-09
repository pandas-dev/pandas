{{ header }}

.. _api.resampling:

==========
Resampling
==========
.. currentmodule:: pandas.api.typing

.. note::
   For an overview, see :ref:`timeseries`.

:class:`pandas.api.typing.Resampler` instances are returned by
resample calls: :func:`pandas.DataFrame.resample`, :func:`pandas.Series.resample`.
The returned object is a :class:`pandas.api.typing.DatetimeIndexResampler`,
:class:`pandas.api.typing.PeriodIndexResampler` or
:class:`pandas.api.typing.TimedeltaIndexResampler`, depending on the type of
the index.

.. autosummary::
   :toctree: api/
   :template: autosummary/class_without_autosummary.rst

   DatetimeIndexResampler
   PeriodIndexResampler
   TimedeltaIndexResampler

Indexing, iteration
~~~~~~~~~~~~~~~~~~~
.. autosummary::
   :toctree: api/

   Resampler.__iter__
   Resampler.groups
   Resampler.indices
   Resampler.get_group

Function application
~~~~~~~~~~~~~~~~~~~~
.. autosummary::
   :toctree: api/

   Resampler.apply
   Resampler.aggregate
   Resampler.transform
   Resampler.pipe

Upsampling
~~~~~~~~~~
.. autosummary::
   :toctree: api/

   Resampler.ffill
   Resampler.bfill
   Resampler.nearest
   Resampler.asfreq
   Resampler.interpolate

Computations / descriptive stats
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
.. autosummary::
   :toctree: api/

   Resampler.count
   Resampler.nunique
   Resampler.first
   Resampler.last
   Resampler.max
   Resampler.mean
   Resampler.median
   Resampler.min
   Resampler.ohlc
   Resampler.prod
   Resampler.size
   Resampler.sem
   Resampler.std
   Resampler.sum
   Resampler.var
   Resampler.quantile
