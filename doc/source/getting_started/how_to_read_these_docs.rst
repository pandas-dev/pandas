.. _howtoread:

{{ header }}

**********************
How to read these docs
**********************

In these tutorials you will see input code inside code blocks such as:

::

    import pandas as pd
    pd.DataFrame({'A': [1, 2, 3]})


or:

.. ipython:: python

    import pandas as pd
    pd.DataFrame({'A': [1, 2, 3]})



The first block is a standard python input, while in the second the ``In [1]:`` indicates the input is inside a `notebook <https://jupyter.org>`__. In Jupyter Notebook the last line is printed and plots are shown inline.

For example:

.. ipython:: python

    a = 1
    a

is equivalent to:

::

    a = 1
    print(a)
