The following tables will be used in the merge examples:

.. ipython:: python

   df1 = pd.DataFrame({"key": ["A", "B", "C", "D"], "value": np.random.randn(4)})
   df1
   df2 = pd.DataFrame({"key": ["B", "D", "D", "E"], "value": np.random.randn(4)})
   df2
