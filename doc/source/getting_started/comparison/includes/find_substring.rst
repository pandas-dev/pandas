Python determines the position of a character in a string with the
:func:`find` function. ``find`` searches for the first position of the
substring. If the substring is found, the function returns its
position. Keep in mind that Python indexes are zero-based and
the function will return -1 if it fails to find the substring.

.. ipython:: python

   tips["sex"].str.find("ale").head()
