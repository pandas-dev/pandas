pandas provides vectorized operations by specifying the individual ``Series`` in the
``DataFrame``. New columns can be assigned in the same way. The :meth:`DataFrame.drop` method drops
a column from the ``DataFrame``.

.. ipython:: python

   tips["total_bill"] = tips["total_bill"] - 2
   tips["new_bill"] = tips["total_bill"] / 2
   tips

   tips = tips.drop("new_bill", axis=1)
