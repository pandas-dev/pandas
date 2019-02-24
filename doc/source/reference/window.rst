{{ header }}

.. _api.window:

======
Window
======
.. currentmodule:: pandas.core.window

Rolling objects are returned by ``.rolling`` calls: :func:`pandas.DataFrame.rolling`, :func:`pandas.Series.rolling`, etc.
Expanding objects are returned by ``.expanding`` calls: :func:`pandas.DataFrame.expanding`, :func:`pandas.Series.expanding`, etc.
EWM objects are returned by ``.ewm`` calls: :func:`pandas.DataFrame.ewm`, :func:`pandas.Series.ewm`, etc.

Standard moving window functions
--------------------------------
.. autosummary::
   :toctree: api/

   Rolling.count
   Rolling.sum
   Rolling.mean
   Rolling.median
   Rolling.var
   Rolling.std
   Rolling.min
   Rolling.max
   Rolling.corr
   Rolling.cov
   Rolling.skew
   Rolling.kurt
   Rolling.apply
   Rolling.aggregate
   Rolling.quantile
   Window.mean
   Window.sum

.. _api.functions_expanding:

Standard expanding window functions
-----------------------------------
.. autosummary::
   :toctree: api/

   Expanding.count
   Expanding.sum
   Expanding.mean
   Expanding.median
   Expanding.var
   Expanding.std
   Expanding.min
   Expanding.max
   Expanding.corr
   Expanding.cov
   Expanding.skew
   Expanding.kurt
   Expanding.apply
   Expanding.aggregate
   Expanding.quantile

Exponentially-weighted moving window functions
----------------------------------------------
.. autosummary::
   :toctree: api/

   EWM.mean
   EWM.std
   EWM.var
   EWM.corr
   EWM.cov
