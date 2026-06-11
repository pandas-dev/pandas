With pandas you can use ``[]`` notation to extract a substring
from a string by position locations. Keep in mind that Python
indexes are zero-based.

.. ipython:: python

   tips["sex"].str[0:1]
