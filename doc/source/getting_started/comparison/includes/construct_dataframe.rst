A pandas ``DataFrame`` can be constructed in many different ways,
but for a small number of values, it is often convenient to specify it as
a Python dictionary, where the keys are the column names
and the values are the data.

.. ipython:: python

   df = pd.DataFrame({"x": [1, 3, 5], "y": [2, 4, 6]})
   df
