:orphan:

The same operation in pandas can be accomplished using
the ``where`` method from ``numpy``.

.. ipython:: python
   :suppress:

   # ensure tips is defined when scanning with flake8-rst
   if 'tips' not in vars():
       tips = {}

.. ipython:: python

   tips["bucket"] = np.where(tips["total_bill"] < 10, "low", "high")
   tips.head()

.. ipython:: python
   :suppress:

   tips = tips.drop("bucket", axis=1)
