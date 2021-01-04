The same operations are expressed in pandas below. Note that these operations do not happen in
place. To make these changes persist, assign the operation back to a variable.

.. ipython:: python

   # keep
   tips[["sex", "total_bill", "tip"]].head()

   # drop
   tips.drop("sex", axis=1).head()

   # rename
   tips.rename(columns={"total_bill": "total_bill_2"}).head()
