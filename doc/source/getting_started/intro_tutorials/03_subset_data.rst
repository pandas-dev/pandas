.. _10min_tut_03_subset:

{{ header }}

.. ipython:: python

    import pandas as pd

.. raw:: html

    <div class="card gs-data">
        <div class="card-header">
            <div class="gs-data-title">
                Data used for this tutorial:
            </div>
        </div>
        <ul class="list-group list-group-flush">
            <li class="list-group-item">

.. include:: includes/titanic.rst

.. ipython:: python

    titanic = pd.read_csv("data/titanic.csv")
    titanic.head()

.. raw:: html

            </li>
        </ul>
    </div>

How do I select a subset of a ``DataFrame``?
============================================

How do I select specific columns from a ``DataFrame``?
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. image:: ../../_static/schemas/03_subset_columns.svg
   :align: center

.. raw:: html

    <ul class="task-bullet">
        <li>

I’m interested in the age of the Titanic passengers.

.. ipython:: python

    ages = titanic["Age"]
    ages.head()

To select a single column, use square brackets ``[]`` with the column
name of the column of interest.

.. raw:: html

        </li>
    </ul>

Each column in a :class:`DataFrame` is a :class:`Series`. As a single column is
selected, the returned object is a pandas :class:`Series`. We can verify this
by checking the type of the output:

.. ipython:: python

    type(titanic["Age"])

And have a look at the ``shape`` of the output:

.. ipython:: python

    titanic["Age"].shape

:attr:`DataFrame.shape` is an attribute (remember :ref:`tutorial on reading and writing <10min_tut_02_read_write>`, do not use parentheses for attributes) of a
pandas ``Series`` and ``DataFrame`` containing the number of rows and
columns: *(nrows, ncolumns)*. A pandas Series is 1-dimensional and only
the number of rows is returned.

.. raw:: html

    <ul class="task-bullet">
        <li>

I’m interested in the age and sex of the Titanic passengers.

.. ipython:: python

    age_sex = titanic[["Age", "Sex"]]
    age_sex.head()

To select multiple columns, use a list of column names within the
selection brackets ``[]``.

.. raw:: html

        </li>
    </ul>

.. note::
    The inner square brackets define a
    :ref:`Python list <python:tut-morelists>` with column names, whereas
    the outer brackets are used to select the data from a pandas
    ``DataFrame`` as seen in the previous example.

The returned data type is a pandas DataFrame:

.. ipython:: python

    type(titanic[["Age", "Sex"]])

.. ipython:: python

    titanic[["Age", "Sex"]].shape

The selection returned a ``DataFrame`` with 891 rows and 2 columns. Remember, a
``DataFrame`` is 2-dimensional with both a row and column dimension.

.. raw:: html

    <div class="d-flex flex-row gs-torefguide">
        <span class="badge badge-info">To user guide</span>

For basic information on indexing, see the user guide section on :ref:`indexing and selecting data <indexing.basics>`.

.. raw:: html

    </div>

How do I filter specific rows from a ``DataFrame``?
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. image:: ../../_static/schemas/03_subset_rows.svg
   :align: center

.. raw:: html

    <ul class="task-bullet">
        <li>

I’m interested in the passengers older than 35 years.

.. ipython:: python

    above_35 = titanic[titanic["Age"] > 35]
    above_35.head()

To select rows based on a conditional expression, use a condition inside
the selection brackets ``[]``.

.. raw:: html

        </li>
    </ul>

The condition inside the selection
brackets ``titanic["Age"] > 35`` checks for which rows the ``Age``
column has a value larger than 35:

.. ipython:: python

    titanic["Age"] > 35

The output of the conditional expression (``>``, but also ``==``,
``!=``, ``<``, ``<=``,… would work) is actually a pandas ``Series`` of
boolean values (either ``True`` or ``False``) with the same number of
rows as the original ``DataFrame``. Such a ``Series`` of boolean values
can be used to filter the ``DataFrame`` by putting it in between the
selection brackets ``[]``. Only rows for which the value is ``True``
will be selected.

We know from before that the original Titanic ``DataFrame`` consists of
891 rows. Let’s have a look at the number of rows which satisfy the
condition by checking the ``shape`` attribute of the resulting
``DataFrame`` ``above_35``:

.. ipython:: python

    above_35.shape

.. raw:: html

    <ul class="task-bullet">
        <li>

I’m interested in the Titanic passengers from cabin class 2 and 3.

.. ipython:: python

    class_23 = titanic[titanic["Pclass"].isin([2, 3])]
    class_23.head()

Similar to the conditional expression, the :func:`~Series.isin` conditional function
returns a ``True`` for each row the values are in the provided list. To
filter the rows based on such a function, use the conditional function
inside the selection brackets ``[]``. In this case, the condition inside
the selection brackets ``titanic["Pclass"].isin([2, 3])`` checks for
which rows the ``Pclass`` column is either 2 or 3.

.. raw:: html

        </li>
    </ul>

The above is equivalent to filtering by rows for which the class is
either 2 or 3 and combining the two statements with an ``|`` (or)
operator:

.. ipython:: python

    class_23 = titanic[(titanic["Pclass"] == 2) | (titanic["Pclass"] == 3)]
    class_23.head()

.. note::
    When combining multiple conditional statements, each condition
    must be surrounded by parentheses ``()``. Moreover, you can not use
    ``or``/``and`` but need to use the ``or`` operator ``|`` and the ``and``
    operator ``&``.

.. raw:: html

    <div class="d-flex flex-row gs-torefguide">
        <span class="badge badge-info">To user guide</span>

See the dedicated section in the user guide about :ref:`boolean indexing <indexing.boolean>` or about the :ref:`isin function <indexing.basics.indexing_isin>`.

.. raw:: html

    </div>

.. raw:: html

    <ul class="task-bullet">
        <li>

I want to work with passenger data for which the age is known.

.. ipython:: python

    age_no_na = titanic[titanic["Age"].notna()]
    age_no_na.head()

The :meth:`~Series.notna` conditional function returns a ``True`` for each row the
values are not an ``Null`` value. As such, this can be combined with the
selection brackets ``[]`` to filter the data table.

.. raw:: html

        </li>
    </ul>

You might wonder what actually changed, as the first 5 lines are still
the same values. One way to verify is to check if the shape has changed:

.. ipython:: python

    age_no_na.shape

.. raw:: html

    <div class="d-flex flex-row gs-torefguide">
        <span class="badge badge-info">To user guide</span>

For more dedicated functions on missing values, see the user guide section about :ref:`handling missing data <missing_data>`.

.. raw:: html

    </div>

.. _10min_tut_03_subset.rows_and_columns:

How do I select specific rows and columns from a ``DataFrame``?
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. image:: ../../_static/schemas/03_subset_columns_rows.svg
   :align: center

.. raw:: html

    <ul class="task-bullet">
        <li>

I’m interested in the names of the passengers older than 35 years.

.. ipython:: python

    adult_names = titanic.loc[titanic["Age"] > 35, "Name"]
    adult_names.head()

In this case, a subset of both rows and columns is made in one go and
just using selection brackets ``[]`` is not sufficient anymore. The
``loc``/``iloc`` operators are required in front of the selection
brackets ``[]``. When using ``loc``/``iloc``, the part before the comma
is the rows you want, and the part after the comma is the columns you
want to select.

.. raw:: html

        </li>
    </ul>

When using the column names, row labels or a condition expression, use
the ``loc`` operator in front of the selection brackets ``[]``. For both
the part before and after the comma, you can use a single label, a list
of labels, a slice of labels, a conditional expression or a colon. Using
a colon specifies you want to select all rows or columns.

.. raw:: html

    <ul class="task-bullet">
        <li>

I’m interested in rows 10 till 25 and columns 3 to 5.

.. ipython:: python

    titanic.iloc[9:25, 2:5]

Again, a subset of both rows and columns is made in one go and just
using selection brackets ``[]`` is not sufficient anymore. When
specifically interested in certain rows and/or columns based on their
position in the table, use the ``iloc`` operator in front of the
selection brackets ``[]``.

.. raw:: html

        </li>
    </ul>

When selecting specific rows and/or columns with ``loc`` or ``iloc``,
new values can be assigned to the selected data. For example, to assign
the name ``anonymous`` to the first 3 elements of the third column:

.. ipython:: python

    titanic.iloc[0:3, 3] = "anonymous"
    titanic.head()

.. raw:: html

    <div class="d-flex flex-row gs-torefguide">
        <span class="badge badge-info">To user guide</span>

See the user guide section on :ref:`different choices for indexing <indexing.choice>` to get more insight in the usage of ``loc`` and ``iloc``.

.. raw:: html

    </div>

.. raw:: html

    <div class="shadow gs-callout gs-callout-remember">
        <h4>REMEMBER</h4>

-  When selecting subsets of data, square brackets ``[]`` are used.
-  Inside these brackets, you can use a single column/row label, a list
   of column/row labels, a slice of labels, a conditional expression or
   a colon.
-  Select specific rows and/or columns using ``loc`` when using the row
   and column names
-  Select specific rows and/or columns using ``iloc`` when using the
   positions in the table
-  You can assign new values to a selection based on ``loc``/``iloc``.

.. raw:: html

    </div>

.. raw:: html

    <div class="d-flex flex-row gs-torefguide">
        <span class="badge badge-info">To user guide</span>

A full overview of indexing is provided in the user guide pages on :ref:`indexing and selecting data <indexing>`.

.. raw:: html

    </div>
