By default, pandas will truncate output of large ``DataFrame``\s to show the first and last rows.
This can be overridden by :ref:`changing the pandas options <options>`, or using
:meth:`DataFrame.head` or :meth:`DataFrame.tail`.

.. ipython:: python

   tips.head(5)
