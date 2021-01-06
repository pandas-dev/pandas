Most pandas operations return copies of the ``Series``/``DataFrame``. To make the changes "stick",
you'll need to either:

* Assign back to the same variable

   .. code-block:: python

      df = df.sort_values("col1")

* Assign to a new variable

   .. code-block:: python

      sorted_df = df.sort_values("col1")

* Use the ``inplace=True`` keyword argument, where available

   .. code-block:: python

      df.sort_values("col1", inplace=True)
