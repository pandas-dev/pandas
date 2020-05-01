.. _10min_tut_07_reshape:

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
                <div data-toggle="collapse" href="#collapsedata" role="button" aria-expanded="false" aria-controls="collapsedata">
                    <span class="badge badge-dark">Titanic data</span>
                </div>
                <div class="collapse" id="collapsedata">
                    <div class="card-body">
                        <p class="card-text">

This tutorial uses the titanic data set, stored as CSV. The data
consists of the following data columns:

-  PassengerId: Id of every passenger.
-  Survived: This feature have value 0 and 1. 0 for not survived and 1
   for survived.
-  Pclass: There are 3 classes: Class 1, Class 2 and Class 3.
-  Name: Name of passenger.
-  Sex: Gender of passenger.
-  Age: Age of passenger.
-  SibSp: Indication that passenger have siblings and spouse.
-  Parch: Whether a passenger is alone or have family.
-  Ticket: Ticket number of passenger.
-  Fare: Indicating the fare.
-  Cabin: The cabin of passenger.
-  Embarked: The embarked category.

.. raw:: html

                        </p>
                    <a href="https://github.com/pandas-dev/pandas/tree/master/doc/data/titanic.csv" class="btn btn-dark btn-sm">To raw data</a>
                </div>
            </div>

.. ipython:: python

    titanic = pd.read_csv("data/titanic.csv")
    titanic.head()

.. raw:: html

        </li>
        <li class="list-group-item">
            <div data-toggle="collapse" href="#collapsedata2" role="button" aria-expanded="false" aria-controls="collapsedata2">
                <span class="badge badge-dark">Air quality data</span>
            </div>
            <div class="collapse" id="collapsedata2">
                <div class="card-body">
                    <p class="card-text">

This tutorial uses air quality data about :math:`NO_2` and Particulate matter less than 2.5
micrometers, made available by
`openaq <https://openaq.org>`__ and using the
`py-openaq <http://dhhagan.github.io/py-openaq/index.html>`__ package.
The ``air_quality_long.csv`` data set provides :math:`NO_2` and
:math:`PM_{25}` values for the measurement stations *FR04014*, *BETR801*
and *London Westminster* in respectively Paris, Antwerp and London.

The air-quality data set has the following columns:

-  city: city where the sensor is used, either Paris, Antwerp or London
-  country: country where the sensor is used, either FR, BE or GB
-  location: the id of the sensor, either *FR04014*, *BETR801* or
   *London Westminster*
-  parameter: the parameter measured by the sensor, either :math:`NO_2`
   or Particulate matter
-  value: the measured value
-  unit: the unit of the measured parameter, in this case ‘µg/m³’

and the index of the ``DataFrame`` is ``datetime``, the datetime of the
measurement.

.. note::
    The air-quality data is provided in a so-called *long format*
    data representation with each observation on a separate row and each
    variable a separate column of the data table. The long/narrow format is
    also known as the `tidy data
    format <https://www.jstatsoft.org/article/view/v059i10>`__.

.. raw:: html

                    </p>
                <a href="https://github.com/pandas-dev/pandas/tree/master/doc/data/air_quality_long.csv" class="btn btn-dark btn-sm">To raw data</a>
            </div>
        </div>

.. ipython:: python

    air_quality = pd.read_csv("data/air_quality_long.csv",
                              index_col="date.utc", parse_dates=True)
    air_quality.head()

.. raw:: html

        </li>
    </ul>
    </div>

How to reshape the layout of tables?
------------------------------------

Sort table rows
~~~~~~~~~~~~~~~

.. raw:: html

    <ul class="task-bullet">
        <li>

I want to sort the titanic data according to the age of the passengers.

.. ipython:: python

    titanic.sort_values(by="Age").head()

.. raw:: html

        </li>
    </ul>

.. raw:: html

    <ul class="task-bullet">
        <li>

I want to sort the titanic data according to the cabin class and age in descending order.

.. ipython:: python

    titanic.sort_values(by=['Pclass', 'Age'], ascending=False).head()

With :meth:`Series.sort_values`, the rows in the table are sorted according to the
defined column(s). The index will follow the row order.

.. raw:: html

        </li>
    </ul>

.. raw:: html

    <div class="d-flex flex-row gs-torefguide">
        <span class="badge badge-info">To user guide</span>

More details about sorting of tables is provided in the using guide section on :ref:`sorting data <basics.sorting>`.

.. raw:: html

   </div>

Long to wide table format
~~~~~~~~~~~~~~~~~~~~~~~~~

Let’s use a small subset of the air quality data set. We focus on
:math:`NO_2` data and only use the first two measurements of each
location (i.e. the head of each group). The subset of data will be
called ``no2_subset``

.. ipython:: python

    # filter for no2 data only
    no2 = air_quality[air_quality["parameter"] == "no2"]

.. ipython:: python

    # use 2 measurements (head) for each location (groupby)
    no2_subset = no2.sort_index().groupby(["location"]).head(2)
    no2_subset

.. image:: ../../_static/schemas/07_pivot.svg
   :align: center

.. raw:: html

    <ul class="task-bullet">
        <li>

I want the values for the three stations as separate columns next to each other

.. ipython:: python

    no2_subset.pivot(columns="location", values="value")

The :meth:`~pandas.pivot_table` function is purely reshaping of the data: a single value
for each index/column combination is required.

.. raw:: html

        </li>
    </ul>

As pandas support plotting of multiple columns (see :ref:`plotting tutorial <10min_tut_04_plotting>`) out of the box, the conversion from
*long* to *wide* table format enables the plotting of the different time
series at the same time:

.. ipython:: python

    no2.head()

.. ipython:: python

    @savefig 7_reshape_columns.png
    no2.pivot(columns="location", values="value").plot()

.. note::
    When the ``index`` parameter is not defined, the existing
    index (row labels) is used.

.. raw:: html

    <div class="d-flex flex-row gs-torefguide">
        <span class="badge badge-info">To user guide</span>

For more information about :meth:`~DataFrame.pivot`, see the user guide section on :ref:`pivoting DataFrame objects <reshaping.reshaping>`.

.. raw:: html

   </div>

Pivot table
~~~~~~~~~~~

.. image:: ../../_static/schemas/07_pivot_table.svg
   :align: center

.. raw:: html

    <ul class="task-bullet">
        <li>

I want the mean concentrations for :math:`NO_2` and :math:`PM_{2.5}` in each of the stations in table form

.. ipython:: python

    air_quality.pivot_table(values="value", index="location",
                            columns="parameter", aggfunc="mean")

In the case of :meth:`~DataFrame.pivot`, the data is only rearranged. When multiple
values need to be aggregated (in this specific case, the values on
different time steps) :meth:`~DataFrame.pivot_table` can be used, providing an
aggregation function (e.g. mean) on how to combine these values.

.. raw:: html

        </li>
    </ul>

Pivot table is a well known concept in spreadsheet software. When
interested in summary columns for each variable separately as well, put
the ``margin`` parameter to ``True``:

.. ipython:: python

    air_quality.pivot_table(values="value", index="location",
                            columns="parameter", aggfunc="mean",
                            margins=True)

.. raw:: html

    <div class="d-flex flex-row gs-torefguide">
        <span class="badge badge-info">To user guide</span>

For more information about :meth:`~DataFrame.pivot_table`, see the user guide section on :ref:`pivot tables <reshaping.pivot>`.

.. raw:: html

   </div>

.. note::
    If case you are wondering, :meth:`~DataFrame.pivot_table` is indeed directly linked
    to :meth:`~DataFrame.groupby`. The same result can be derived by grouping on both
    ``parameter`` and ``location``:

    ::

        air_quality.groupby(["parameter", "location"]).mean()

.. raw:: html

    <div class="d-flex flex-row gs-torefguide">
        <span class="badge badge-info">To user guide</span>

Have a look at :meth:`~DataFrame.groupby` in combination with :meth:`~DataFrame.unstack` at the user guide section on :ref:`combining stats and groupby <reshaping.combine_with_groupby>`.

.. raw:: html

   </div>

Wide to long format
~~~~~~~~~~~~~~~~~~~

Starting again from the wide format table created in the previous
section:

.. ipython:: python

    no2_pivoted = no2.pivot(columns="location", values="value").reset_index()
    no2_pivoted.head()

.. image:: ../../_static/schemas/07_melt.svg
   :align: center

.. raw:: html

    <ul class="task-bullet">
        <li>

I want to collect all air quality :math:`NO_2` measurements in a single column (long format)

.. ipython:: python

    no_2 = no2_pivoted.melt(id_vars="date.utc")
    no_2.head()

The :func:`pandas.melt` method on a ``DataFrame`` converts the data table from wide
format to long format. The column headers become the variable names in a
newly created column.

.. raw:: html

        </li>
    </ul>

The solution is the short version on how to apply :func:`pandas.melt`. The method
will *melt* all columns NOT mentioned in ``id_vars`` together into two
columns: A columns with the column header names and a column with the
values itself. The latter column gets by default the name ``value``.

The :func:`pandas.melt` method can be defined in more detail:

.. ipython:: python

    no_2 = no2_pivoted.melt(id_vars="date.utc",
                            value_vars=["BETR801",
                                        "FR04014",
                                        "London Westminster"],
                            value_name="NO_2",
                            var_name="id_location")
    no_2.head()

The result in the same, but in more detail defined:

-  ``value_vars`` defines explicitly which columns to *melt* together
-  ``value_name`` provides a custom column name for the values column
   instead of the default columns name ``value``
-  ``var_name`` provides a custom column name for the columns collecting
   the column header names. Otherwise it takes the index name or a
   default ``variable``

Hence, the arguments ``value_name`` and ``var_name`` are just
user-defined names for the two generated columns. The columns to melt
are defined by ``id_vars`` and ``value_vars``.

.. raw:: html

    <div class="d-flex flex-row gs-torefguide">
        <span class="badge badge-info">To user guide</span>

Conversion from wide to long format with :func:`pandas.melt` is explained in the user guide section on :ref:`reshaping by melt <reshaping.melt>`.

.. raw:: html

   </div>

.. raw:: html

    <div class="shadow gs-callout gs-callout-remember">
        <h4>REMEMBER</h4>

-  Sorting by one or more columns is supported by ``sort_values``
-  The ``pivot`` function is purely restructering of the data,
   ``pivot_table`` supports aggregations
-  The reverse of ``pivot`` (long to wide format) is ``melt`` (wide to
   long format)

.. raw:: html

   </div>

.. raw:: html

    <div class="d-flex flex-row gs-torefguide">
        <span class="badge badge-info">To user guide</span>

A full overview is available in the user guide on the pages about :ref:`reshaping and pivoting <reshaping>`.

.. raw:: html

   </div>
