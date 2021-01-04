Python determines the length of a character string with the ``len`` function.
In Python 3, all strings are Unicode strings. ``len`` includes trailing blanks.
Use ``len`` and ``rstrip`` to exclude trailing blanks.

.. ipython:: python

   tips["time"].str.len().head()
   tips["time"].str.rstrip().str.len().head()
