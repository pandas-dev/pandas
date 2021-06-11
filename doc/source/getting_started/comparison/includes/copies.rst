Most pandas operations return copies of the ``Series``/``DataFrame``. To make the changes "stick",
you'll need to either assign to a new variable:

   .. code-block:: python

      sorted_df = df.sort_values("col1")


or overwrite the original one:

   .. code-block:: python

      df = df.sort_values("col1")

.. note::

   You will see an ``inplace=True`` keyword argument available for some methods:

   .. code-block:: python

      df.sort_values("col1", inplace=True)

   Its use is discouraged. :ref:`More information. <indexing.view_versus_copy>`
