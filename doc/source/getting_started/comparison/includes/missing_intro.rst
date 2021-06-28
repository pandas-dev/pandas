pandas represents missing data with the special float value ``NaN`` (not a number).  Many of the
semantics are the same; for example missing data propagates through numeric operations, and is
ignored by default for aggregations.

.. ipython:: python

   outer_join
   outer_join["value_x"] + outer_join["value_y"]
   outer_join["value_x"].sum()
