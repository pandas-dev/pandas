.. _future_udf_behavior:
:orphan:

{{ header }}

*******************
Future UDF Behavior
*******************

pandas is experimenting with improving the behavior of methods that take a
user-defined function (UDF). These methods include ``.apply``, ``.agg``, ``.transform``,
and ``.filter``. The goal is to make these methods behave in a more predictable
and consistent manner, reducing the complexity of their implementation, and improving
performance where possible. This page details the differences between the old and
new behaviors, as well as providing some context behind each change that is being made.

There are a great number of changes that are planned. In order to transition in a
reasonable manner for users, all changes are behind an experimental "future_udf_behavior"
option. This is currently experimental and subject to breaking changes without notice.
Users can opt into the new behavior and provide feedback. Once the improvements have
been made, this option will be declared no longer experimental. pandas will then raise
a ``FutureWarning`` that the default value of this option will be set to ``True`` in
a future version. Once the default is ``True``, users can still override it to ``False``.
After a sufficient amount of time, pandas will remove this option altogether and only
the future behavior will remain.

``DataFrame.agg`` with list-likes
---------------------------------

Previously, using ``DataFrame.agg`` with a list-like argument would transpose the result when
compared with just providing a single aggregation function.

.. ipython:: python

   df = pd.DataFrame({'a': [1, 2, 3], 'b': [4, 5, 6], 'c': [7, 8, 9]})

   df.agg('sum')
   df.agg(['sum'])

This transpose no longer occurs, making the result more consistent.

.. ipython:: python

   with pd.option_context('future_udf_behavior', True):
       result = df.agg(['sum'])
   result

   with pd.option_context('future_udf_behavior', True):
       result = df.agg(['sum', 'mean'])
   result

``DataFrame.groupby(...).agg`` with list-likes
----------------------------------------------

Previously, using ``DataFrame.groupby(...).agg`` with a list-like argument would put the
columns as the first level of the resulting hierarchical columns. The result is
that the columns for each aggregation function are separated, inconsistent with the result
for a single aggregator.

.. ipython:: python

   df.groupby("a").agg('sum')
   df.groupby("a").agg(["sum", "min"])

Now the levels are swapped, so that the columns for each aggregation are together.

.. ipython:: python

   with pd.option_context('future_udf_behavior', True):
       result = df.groupby("a").agg(["sum", "min"])
   result
