:orphan:

The same operation in pandas can be accomplished using
the ``where`` method from ``numpy``.

.. ipython:: python
   :flake8-add-ignore: F821

   tips["bucket"] = np.where(tips["total_bill"] < 10, "low", "high")
   tips.head()

.. ipython:: python
   :flake8-add-ignore: F821
   :suppress:

   tips = tips.drop("bucket", axis=1)
