Most pandas operations return copies of the ``Series``/``DataFrame``. To make the changes "stick",
you'll need to either assign to a new variable:

   .. code-block:: python

      sorted_df = df.sort_values("col1")


or overwrite the original one:

   .. code-block:: python

      df = df.sort_values("col1")

.. note::

   You will see an ``inplace=True`` or ``copy=False`` keyword argument available for
   some methods:

   .. code-block:: python

      df.replace(5, inplace=True)

   There is an active discussion about deprecating and removing ``inplace`` and ``copy`` for
   most methods (e.g. ``dropna``) except for a very small subset of methods
   (including ``replace``). Both keywords won't be
   necessary anymore in the context of Copy-on-Write. The proposal can be found
   `here <https://github.com/pandas-dev/pandas/pull/51466>`_.
