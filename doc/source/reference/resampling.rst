{{ header }}

.. _api.resampling:

==========
Resampling
==========
.. currentmodule:: pandas.core.resample

Resampler objects are returned by resample calls: :func:`pandas.DataFrame.resample`, :func:`pandas.Series.resample`.

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
   Resampler.backfill
   Resampler.bfill
   Resampler.pad
   Resampler.nearest
   Resampler.fillna
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
