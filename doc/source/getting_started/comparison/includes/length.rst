You can find the length of a character string with :meth:`Series.str.len`.
In Python 3, all strings are Unicode strings. ``len`` includes trailing blanks.
Use ``len`` and ``rstrip`` to exclude trailing blanks.

.. ipython:: python

   tips["time"].str.len()
   tips["time"].str.rstrip().str.len()
