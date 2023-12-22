The same operations are expressed in pandas below.

Keep certain columns
''''''''''''''''''''

.. ipython:: python

   tips[["sex", "total_bill", "tip"]]

Drop a column
'''''''''''''

.. ipython:: python

   tips.drop("sex", axis=1)

Rename a column
'''''''''''''''

.. ipython:: python

   tips.rename(columns={"total_bill": "total_bill_2"})
