:orphan:

pandas objects have a :meth:`DataFrame.sort_values` method, which
takes a list of columns to sort by.

.. ipython:: python
   :flake8-add-ignore: F821

   tips = tips.sort_values(["sex", "total_bill"])
   tips.head()
