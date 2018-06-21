.. _release:

.. currentmodule:: pandas

.. ipython:: python
   :suppress:

   import pandas as pd
   import numpy as np
   np.random.seed(123456)
   np.set_printoptions(precision=4, suppress=True)
   import matplotlib.pyplot as plt
   plt.close('all')

   pd.options.display.max_rows=15
   import pandas.util.testing as tm

*************
Release Notes
*************

This is the list of changes to pandas between each release. For full details,
see the commit logs at http://github.com/pandas-dev/pandas. For install and
upgrade instructions, see :ref:`install`.

0.23 release
------------

.. toctree::
   :maxdepth: 1

   /source/whatsnew/v0.23.2.txt
   whatsnew/v0.23.1.txt
   whatsnew/v0.23.0.txt
