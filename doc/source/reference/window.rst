{{ header }}

.. _api.window:

======
Window
======

Rolling objects are returned by ``.rolling`` calls: :func:`pandas.DataFrame.rolling`, :func:`pandas.Series.rolling`, etc.
Expanding objects are returned by ``.expanding`` calls: :func:`pandas.DataFrame.expanding`, :func:`pandas.Series.expanding`, etc.
ExponentialMovingWindow objects are returned by ``.ewm`` calls: :func:`pandas.DataFrame.ewm`, :func:`pandas.Series.ewm`, etc.

.. _api.functions_rolling:

Rolling window functions
------------------------
.. currentmodule:: pandas.core.window.rolling

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
   Rolling.sem

.. _api.functions_window:

Weighted window functions
-------------------------
.. currentmodule:: pandas.core.window.rolling

.. autosummary::
   :toctree: api/

   Window.mean
   Window.sum
   Window.var
   Window.std

.. _api.functions_expanding:

Expanding window functions
--------------------------
.. currentmodule:: pandas.core.window.expanding

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
   Expanding.sem

.. _api.functions_ewm:

Exponentially-weighted window functions
---------------------------------------
.. currentmodule:: pandas.core.window.ewm

.. autosummary::
   :toctree: api/

   ExponentialMovingWindow.mean
   ExponentialMovingWindow.std
   ExponentialMovingWindow.var
   ExponentialMovingWindow.corr
   ExponentialMovingWindow.cov

.. _api.indexers_window:

Window indexer
--------------
.. currentmodule:: pandas

Base class for defining custom window boundaries.

.. autosummary::
   :toctree: api/

   api.indexers.BaseIndexer
   api.indexers.FixedForwardWindowIndexer
   api.indexers.VariableOffsetWindowIndexer
