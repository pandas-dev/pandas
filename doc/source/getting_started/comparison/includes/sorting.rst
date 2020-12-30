:orphan:

pandas objects have a :meth:`DataFrame.sort_values` method, which
takes a list of columns to sort by.

.. ipython:: python
   :suppress:

   # ensure tips is defined when scanning with flake8-rst
   if 'tips' not in vars():
       tips = {}

.. ipython:: python

   tips = tips.sort_values(["sex", "total_bill"])
   tips.head()
