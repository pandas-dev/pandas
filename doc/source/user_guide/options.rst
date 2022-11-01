.. _options:

{{ header }}

********************
Options and settings
********************

Overview
--------
pandas has an options API configure and customize global behavior related to
:class:`DataFrame` display, data behavior and more.

Options have a full "dotted-style", case-insensitive name (e.g. ``display.max_rows``).
You can get/set options directly as attributes of the top-level ``options`` attribute:

.. ipython:: python

   import pandas as pd

   pd.options.display.max_rows
   pd.options.display.max_rows = 999
   pd.options.display.max_rows

The API is composed of 5 relevant functions, available directly from the ``pandas``
namespace:

* :func:`~pandas.get_option` / :func:`~pandas.set_option` - get/set the value of a single option.
* :func:`~pandas.reset_option` - reset one or more options to their default value.
* :func:`~pandas.describe_option` - print the descriptions of one or more options.
* :func:`~pandas.option_context` - execute a codeblock with a set of options
  that revert to prior settings after execution.

.. note::

   Developers can check out `pandas/core/config_init.py <https://github.com/pandas-dev/pandas/blob/main/pandas/core/config_init.py>`_ for more information.

All of the functions above accept a regexp pattern (``re.search`` style) as an argument,
to match an unambiguous substring:

.. ipython:: python

   pd.get_option("display.chop_threshold")
   pd.set_option("display.chop_threshold", 2)
   pd.get_option("display.chop_threshold")
   pd.set_option("chop", 4)
   pd.get_option("display.chop_threshold")


The following will **not work** because it matches multiple option names, e.g.
``display.max_colwidth``, ``display.max_rows``, ``display.max_columns``:

.. ipython:: python
   :okexcept:

   pd.get_option("max")


.. warning::

    Using this form of shorthand may cause your code to break if new options with similar names are added in future versions.


.. ipython:: python
   :suppress:
   :okwarning:

   pd.reset_option("all")

.. _options.available:

Available options
-----------------

You can get a list of available options and their descriptions with :func:`~pandas.describe_option`. When called
with no argument :func:`~pandas.describe_option` will print out the descriptions for all available options.

.. ipython:: python

   pd.describe_option()

Getting and setting options
---------------------------

As described above, :func:`~pandas.get_option` and :func:`~pandas.set_option`
are available from the pandas namespace.  To change an option, call
``set_option('option regex', new_value)``.

.. ipython:: python

   pd.get_option("mode.sim_interactive")
   pd.set_option("mode.sim_interactive", True)
   pd.get_option("mode.sim_interactive")

.. note::

   The option ``'mode.sim_interactive'`` is mostly used for debugging purposes.

You can use :func:`~pandas.reset_option` to revert to a setting's default value

.. ipython:: python
   :suppress:

   pd.reset_option("display.max_rows")

.. ipython:: python

   pd.get_option("display.max_rows")
   pd.set_option("display.max_rows", 999)
   pd.get_option("display.max_rows")
   pd.reset_option("display.max_rows")
   pd.get_option("display.max_rows")


It's also possible to reset multiple options at once (using a regex):

.. ipython:: python
   :okwarning:

   pd.reset_option("^display")


:func:`~pandas.option_context` context manager has been exposed through
the top-level API, allowing you to execute code with given option values. Option values
are restored automatically when you exit the ``with`` block:

.. ipython:: python

   with pd.option_context("display.max_rows", 10, "display.max_columns", 5):
       print(pd.get_option("display.max_rows"))
       print(pd.get_option("display.max_columns"))
   print(pd.get_option("display.max_rows"))
   print(pd.get_option("display.max_columns"))


Setting startup options in Python/IPython environment
-----------------------------------------------------

Using startup scripts for the Python/IPython environment to import pandas and set options makes working with pandas more efficient.
To do this, create a ``.py`` or ``.ipy`` script in the startup directory of the desired profile.
An example where the startup folder is in a default IPython profile can be found at:

.. code-block:: none

  $IPYTHONDIR/profile_default/startup

More information can be found in the `IPython documentation
<https://ipython.org/ipython-doc/stable/interactive/tutorial.html#startup-files>`__.  An example startup script for pandas is displayed below:

.. code-block:: python

  import pandas as pd

  pd.set_option("display.max_rows", 999)
  pd.set_option("display.precision", 5)

.. _options.frequently_used:

Frequently used options
-----------------------
The following is a demonstrates the more frequently used display options.

``display.max_rows`` and ``display.max_columns`` sets the maximum number
of rows and columns displayed when a frame is pretty-printed. Truncated
lines are replaced by an ellipsis.

.. ipython:: python

   df = pd.DataFrame(np.random.randn(7, 2))
   pd.set_option("display.max_rows", 7)
   df
   pd.set_option("display.max_rows", 5)
   df
   pd.reset_option("display.max_rows")

Once the ``display.max_rows`` is exceeded, the ``display.min_rows`` options
determines how many rows are shown in the truncated repr.

.. ipython:: python

   pd.set_option("display.max_rows", 8)
   pd.set_option("display.min_rows", 4)
   # below max_rows -> all rows shown
   df = pd.DataFrame(np.random.randn(7, 2))
   df
   # above max_rows -> only min_rows (4) rows shown
   df = pd.DataFrame(np.random.randn(9, 2))
   df
   pd.reset_option("display.max_rows")
   pd.reset_option("display.min_rows")

``display.expand_frame_repr`` allows for the representation of a
:class:`DataFrame` to stretch across pages, wrapped over the all the columns.

.. ipython:: python

   df = pd.DataFrame(np.random.randn(5, 10))
   pd.set_option("expand_frame_repr", True)
   df
   pd.set_option("expand_frame_repr", False)
   df
   pd.reset_option("expand_frame_repr")

``display.large_repr`` displays a :class:`DataFrame` that exceed
``max_columns`` or ``max_rows`` as a truncated frame or summary.

.. ipython:: python

   df = pd.DataFrame(np.random.randn(10, 10))
   pd.set_option("display.max_rows", 5)
   pd.set_option("large_repr", "truncate")
   df
   pd.set_option("large_repr", "info")
   df
   pd.reset_option("large_repr")
   pd.reset_option("display.max_rows")

``display.max_colwidth`` sets the maximum width of columns.  Cells
of this length or longer will be truncated with an ellipsis.

.. ipython:: python

   df = pd.DataFrame(
       np.array(
           [
               ["foo", "bar", "bim", "uncomfortably long string"],
               ["horse", "cow", "banana", "apple"],
           ]
       )
   )
   pd.set_option("max_colwidth", 40)
   df
   pd.set_option("max_colwidth", 6)
   df
   pd.reset_option("max_colwidth")

``display.max_info_columns`` sets a threshold for the number of columns
displayed when calling :meth:`~pandas.DataFrame.info`.

.. ipython:: python

   df = pd.DataFrame(np.random.randn(10, 10))
   pd.set_option("max_info_columns", 11)
   df.info()
   pd.set_option("max_info_columns", 5)
   df.info()
   pd.reset_option("max_info_columns")

``display.max_info_rows``: :meth:`~pandas.DataFrame.info` will usually show null-counts for each column.
For a large :class:`DataFrame`, this can be quite slow. ``max_info_rows`` and ``max_info_cols``
limit this null check to the specified rows and columns respectively. The :meth:`~pandas.DataFrame.info`
keyword argument ``show_counts=True`` will override this.

.. ipython:: python

   df = pd.DataFrame(np.random.choice([0, 1, np.nan], size=(10, 10)))
   df
   pd.set_option("max_info_rows", 11)
   df.info()
   pd.set_option("max_info_rows", 5)
   df.info()
   pd.reset_option("max_info_rows")

``display.precision`` sets the output display precision in terms of decimal places.

.. ipython:: python

   df = pd.DataFrame(np.random.randn(5, 5))
   pd.set_option("display.precision", 7)
   df
   pd.set_option("display.precision", 4)
   df

``display.chop_threshold`` sets the rounding threshold to zero when displaying a
:class:`Series` or :class:`DataFrame`. This setting does not change the
precision at which the number is stored.

.. ipython:: python

   df = pd.DataFrame(np.random.randn(6, 6))
   pd.set_option("chop_threshold", 0)
   df
   pd.set_option("chop_threshold", 0.5)
   df
   pd.reset_option("chop_threshold")

``display.colheader_justify`` controls the justification of the headers.
The options are ``'right'``, and ``'left'``.

.. ipython:: python

   df = pd.DataFrame(
       np.array([np.random.randn(6), np.random.randint(1, 9, 6) * 0.1, np.zeros(6)]).T,
       columns=["A", "B", "C"],
       dtype="float",
   )
   pd.set_option("colheader_justify", "right")
   df
   pd.set_option("colheader_justify", "left")
   df
   pd.reset_option("colheader_justify")


.. _basics.console_output:

Number formatting
------------------

pandas also allows you to set how numbers are displayed in the console.
This option is not set through the ``set_options`` API.

Use the ``set_eng_float_format`` function
to alter the floating-point formatting of pandas objects to produce a particular
format.

.. ipython:: python

   import numpy as np

   pd.set_eng_float_format(accuracy=3, use_eng_prefix=True)
   s = pd.Series(np.random.randn(5), index=["a", "b", "c", "d", "e"])
   s / 1.0e3
   s / 1.0e6

.. ipython:: python
   :suppress:
   :okwarning:

   pd.reset_option("^display")

Use :meth:`~pandas.DataFrame.round` to specifically control rounding of an individual :class:`DataFrame`

.. _options.east_asian_width:

Unicode formatting
------------------

.. warning::

   Enabling this option will affect the performance for printing of DataFrame and Series (about 2 times slower).
   Use only when it is actually required.

Some East Asian countries use Unicode characters whose width corresponds to two Latin characters.
If a DataFrame or Series contains these characters, the default output mode may not align them properly.

.. ipython:: python

   df = pd.DataFrame({"国籍": ["UK", "日本"], "名前": ["Alice", "しのぶ"]})
   df

Enabling ``display.unicode.east_asian_width`` allows pandas to check each character's "East Asian Width" property.
These characters can be aligned properly by setting this option to ``True``. However, this will result in longer render
times than the standard ``len`` function.

.. ipython:: python

   pd.set_option("display.unicode.east_asian_width", True)
   df

In addition, Unicode characters whose width is "ambiguous" can either be 1 or 2 characters wide depending on the
terminal setting or encoding. The option ``display.unicode.ambiguous_as_wide`` can be used to handle the ambiguity.

By default, an "ambiguous" character's width, such as "¡" (inverted exclamation) in the example below, is taken to be 1.

.. ipython:: python

   df = pd.DataFrame({"a": ["xxx", "¡¡"], "b": ["yyy", "¡¡"]})
   df


Enabling ``display.unicode.ambiguous_as_wide`` makes pandas interpret these characters' widths to be 2.
(Note that this option will only be effective when ``display.unicode.east_asian_width`` is enabled.)

However, setting this option incorrectly for your terminal will cause these characters to be aligned incorrectly:

.. ipython:: python

   pd.set_option("display.unicode.ambiguous_as_wide", True)
   df


.. ipython:: python
   :suppress:

   pd.set_option("display.unicode.east_asian_width", False)
   pd.set_option("display.unicode.ambiguous_as_wide", False)

.. _options.table_schema:

Table schema display
--------------------

:class:`DataFrame` and :class:`Series` will publish a Table Schema representation
by default. This can be enabled globally with the
``display.html.table_schema`` option:

.. ipython:: python

  pd.set_option("display.html.table_schema", True)

Only ``'display.max_rows'`` are serialized and published.


.. ipython:: python
    :suppress:

    pd.reset_option("display.html.table_schema")
