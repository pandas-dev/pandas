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

   The ``copy`` keyword has been deprecated in most methods as part of
   Copy-on-Write, which is the default behavior in pandas 3.0. The ``inplace``
   keyword is still available for a subset of methods (including ``replace``).
