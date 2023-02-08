Most pandas operations return copies of the ``Series``/``DataFrame``. To make the changes "stick",
you'll need to either assign to a new variable:

   .. code-block:: python

      sorted_df = df.sort_values("col1")


or overwrite the original one:

   .. code-block:: python

      df = df.sort_values("col1")
