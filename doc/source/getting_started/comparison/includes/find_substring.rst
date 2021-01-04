You can find the position of a character in a column of strings with the :meth:`Series.str.find`
method. ``find`` searches for the first position of the substring. If the substring is found, the
method returns its position. If not found, it returns ``-1``. Keep in mind that Python indexes are
zero-based.

.. ipython:: python

   tips["sex"].str.find("ale")
