pandas has a :meth:`DataFrame.sort_values` method, which takes a list of columns to sort by.

.. ipython:: python

   tips = tips.sort_values(["sex", "total_bill"])
   tips
