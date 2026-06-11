pandas provides a flexible ``groupby`` mechanism that allows similar aggregations. See the
:ref:`groupby documentation<groupby>` for more details and examples.

.. ipython:: python

   tips_summed = tips.groupby(["sex", "smoker"])[["total_bill", "tip"]].sum()
   tips_summed
