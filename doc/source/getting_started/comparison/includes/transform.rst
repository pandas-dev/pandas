pandas provides a :ref:`groupby.transform` mechanism that allows these type of operations to be
succinctly expressed in one operation.

.. ipython:: python

   gb = tips.groupby("smoker")["total_bill"]
   tips["adj_total_bill"] = tips["total_bill"] - gb.transform("mean")
   tips
