.. ipython:: python

   tips["date1"] = pd.Timestamp("2013-01-15")
   tips["date2"] = pd.Timestamp("2015-02-15")
   tips["date1_year"] = tips["date1"].dt.year
   tips["date2_month"] = tips["date2"].dt.month
   tips["date1_next"] = tips["date1"] + pd.offsets.MonthBegin()
   tips["months_between"] = tips["date2"].dt.to_period("M") - tips[
       "date1"
   ].dt.to_period("M")

   tips[
       ["date1", "date2", "date1_year", "date2_month", "date1_next", "months_between"]
   ]

.. ipython:: python
   :suppress:

   tips = tips.drop(
       ["date1", "date2", "date1_year", "date2_month", "date1_next", "months_between"],
       axis=1,
   )
