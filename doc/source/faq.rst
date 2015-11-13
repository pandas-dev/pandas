.. currentmodule:: pandas
.. _faq:

********************************
Frequently Asked Questions (FAQ)
********************************

.. ipython:: python
   :suppress:

   import numpy as np
   np.random.seed(123456)
   np.set_printoptions(precision=4, suppress=True)
   import pandas as pd
   pd.options.display.max_rows = 15
   import matplotlib
   try:
      matplotlib.style.use('ggplot')
   except AttributeError:
      pd.options.display.mpl_style = 'default'
   import matplotlib.pyplot as plt
   plt.close('all')

.. _df-memory-usage:

DataFrame memory usage
----------------------
As of pandas version 0.15.0, the memory usage of a dataframe (including
the index) is shown when accessing the ``info`` method of a dataframe. A
configuration option, ``display.memory_usage`` (see :ref:`options`),
specifies if the dataframe's memory usage will be displayed when
invoking the ``df.info()`` method.

For example, the memory usage of the dataframe below is shown
when calling ``df.info()``:

.. ipython:: python

    dtypes = ['int64', 'float64', 'datetime64[ns]', 'timedelta64[ns]',
              'complex128', 'object', 'bool']
    n = 5000
    data = dict([ (t, np.random.randint(100, size=n).astype(t))
                    for t in dtypes])
    df = pd.DataFrame(data)
    df['categorical'] = df['object'].astype('category')

    df.info()

The ``+`` symbol indicates that the true memory usage could be higher, because
pandas does not count the memory used by values in columns with
``dtype=object``.

.. versionadded:: 0.17.1

Passing ``memory_usage='deep'`` will enable a more accurate memory usage report,
that accounts for the full usage of the contained objects. This is optional
as it can be expensive to do this deeper introspection.

.. ipython:: python

   df.info(memory_usage='deep')

By default the display option is set to ``True`` but can be explicitly
overridden by passing the ``memory_usage`` argument when invoking ``df.info()``.

The memory usage of each column can be found by calling the ``memory_usage``
method. This returns a Series with an index represented by column names
and memory usage of each column shown in bytes. For the dataframe above,
the memory usage of each column and the total memory usage of the
dataframe can be found with the memory_usage method:

.. ipython:: python

    df.memory_usage()

    # total memory usage of dataframe
    df.memory_usage().sum()

By default the memory usage of the dataframe's index is not shown in the
returned Series, the memory usage of the index can be shown by passing
the ``index=True`` argument:

.. ipython:: python

    df.memory_usage(index=True)

The memory usage displayed by the ``info`` method utilizes the
``memory_usage`` method to determine the memory usage of a dataframe
while also formatting the output in human-readable units (base-2
representation; i.e., 1KB = 1024 bytes).

See also :ref:`Categorical Memory Usage <categorical.memory>`.

Byte-Ordering Issues
--------------------
Occasionally you may have to deal with data that were created on a machine with
a different byte order than the one on which you are running Python. To deal
with this issue you should convert the underlying NumPy array to the native
system byte order *before* passing it to Series/DataFrame/Panel constructors
using something similar to the following:

.. ipython:: python

   x = np.array(list(range(10)), '>i4') # big endian
   newx = x.byteswap().newbyteorder() # force native byteorder
   s = pd.Series(newx)

See `the NumPy documentation on byte order
<http://docs.scipy.org/doc/numpy/user/basics.byteswapping.html>`__ for more
details.


Visualizing Data in Qt applications
-----------------------------------

.. warning::

    The ``qt`` support is **deprecated and will be removed in a future version**.
    We refer users to the external package `pandas-qt <https://github.com/datalyze-solutions/pandas-qt>`_.

There is experimental support for visualizing DataFrames in PyQt4 and PySide
applications. At the moment you can display and edit the values of the cells
in the DataFrame. Qt will take care of displaying just the portion of the
DataFrame that is currently visible and the edits will be immediately saved to
the underlying DataFrame

To demonstrate this we will create a simple PySide application that will switch
between two editable DataFrames. For this will use the ``DataFrameModel`` class
that handles the access to the DataFrame, and the ``DataFrameWidget``, which is
just a thin layer around the ``QTableView``.

.. code-block:: python

	import numpy as np
	import pandas as pd
	from pandas.sandbox.qtpandas import DataFrameModel, DataFrameWidget
	from PySide import QtGui, QtCore

	# Or if you use PyQt4:
	# from PyQt4 import QtGui, QtCore

	class MainWidget(QtGui.QWidget):
	    def __init__(self, parent=None):
	        super(MainWidget, self).__init__(parent)

	        # Create two DataFrames
	        self.df1 = pd.DataFrame(np.arange(9).reshape(3, 3),
	                                columns=['foo', 'bar', 'baz'])
	        self.df2 = pd.DataFrame({
	                'int': [1, 2, 3],
	                'float': [1.5, 2.5, 3.5],
	                'string': ['a', 'b', 'c'],
	                'nan': [np.nan, np.nan, np.nan]
	            }, index=['AAA', 'BBB', 'CCC'],
	            columns=['int', 'float', 'string', 'nan'])

	        # Create the widget and set the first DataFrame
	        self.widget = DataFrameWidget(self.df1)

	        # Create the buttons for changing DataFrames
	        self.button_first = QtGui.QPushButton('First')
	        self.button_first.clicked.connect(self.on_first_click)
	        self.button_second = QtGui.QPushButton('Second')
	        self.button_second.clicked.connect(self.on_second_click)

	        # Set the layout
	        vbox = QtGui.QVBoxLayout()
	        vbox.addWidget(self.widget)
	        hbox = QtGui.QHBoxLayout()
	        hbox.addWidget(self.button_first)
	        hbox.addWidget(self.button_second)
	        vbox.addLayout(hbox)
	        self.setLayout(vbox)

	    def on_first_click(self):
	    	'''Sets the first DataFrame'''
	        self.widget.setDataFrame(self.df1)

	    def on_second_click(self):
	    	'''Sets the second DataFrame'''
	        self.widget.setDataFrame(self.df2)

	if __name__ == '__main__':
	    import sys

	    # Initialize the application
	    app = QtGui.QApplication(sys.argv)
	    mw = MainWidget()
	    mw.show()
	    app.exec_()
