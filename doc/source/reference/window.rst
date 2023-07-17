{{ header }}

.. _api.window:

======
Window
======

:class:`pandas.api.typing.Rolling` instances are returned by ``.rolling`` calls:
:func:`pandas.DataFrame.rolling` and :func:`pandas.Series.rolling`.
:class:`pandas.api.typing.Expanding` instances are returned by ``.expanding`` calls:
:func:`pandas.DataFrame.expanding` and :func:`pandas.Series.expanding`.
:class:`pandas.api.typing.ExponentialMovingWindow` instances are returned by ``.ewm``
calls: :func:`pandas.DataFrame.ewm` and :func:`pandas.Series.ewm`.

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
   Rolling.rank

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
   Expanding.rank

.. _api.functions_ewm:

Exponentially-weighted window functions
---------------------------------------
.. currentmodule:: pandas.core.window.ewm

.. autosummary::
   :toctree: api/

   ExponentialMovingWindow.mean
   ExponentialMovingWindow.sum
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
