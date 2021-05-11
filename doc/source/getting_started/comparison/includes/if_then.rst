The same operation in pandas can be accomplished using
the ``where`` method from ``numpy``.

.. ipython:: python

   tips["bucket"] = np.where(tips["total_bill"] < 10, "low", "high")
   tips

.. ipython:: python
   :suppress:

   tips = tips.drop("bucket", axis=1)
