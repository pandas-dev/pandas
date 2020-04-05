.. _10min_tut_08_combine:

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
                    <span class="badge badge-dark">Air quality Nitrate data</span>
                </div>
                <div class="collapse" id="collapsedata">
                    <div class="card-body">
                        <p class="card-text">

For this tutorial, air quality data about :math:`NO_2` is used, made available by
`openaq <https://openaq.org>`__ and downloaded using the
`py-openaq <http://dhhagan.github.io/py-openaq/index.html>`__ package.

The ``air_quality_no2_long.csv`` data set provides :math:`NO_2`
values for the measurement stations *FR04014*, *BETR801* and *London
Westminster* in respectively Paris, Antwerp and London.

.. raw:: html

                        </p>
                    <a href="https://github.com/pandas-dev/pandas/tree/master/doc/data/air_quality_no2_long.csv" class="btn btn-dark btn-sm">To raw data</a>
                </div>
            </div>

.. ipython:: python

    air_quality_no2 = pd.read_csv("data/air_quality_no2_long.csv",
                                  parse_dates=True)
    air_quality_no2 = air_quality_no2[["date.utc", "location",
                                       "parameter", "value"]]
    air_quality_no2.head()

.. raw:: html

        </li>
        <li class="list-group-item">
            <div data-toggle="collapse" href="#collapsedata2" role="button" aria-expanded="false" aria-controls="collapsedata2">
                <span class="badge badge-dark">Air quality Particulate matter data</span>
            </div>
            <div class="collapse" id="collapsedata2">
                <div class="card-body">
                    <p class="card-text">

For this tutorial, air quality data about Particulate
matter less than 2.5 micrometers is used, made available by
`openaq <https://openaq.org>`__ and downloaded using the
`py-openaq <http://dhhagan.github.io/py-openaq/index.html>`__ package.

The ``air_quality_pm25_long.csv`` data set provides :math:`PM_{25}`
values for the measurement stations *FR04014*, *BETR801* and *London
Westminster* in respectively Paris, Antwerp and London.

.. raw:: html

                    </p>
                <a href="https://github.com/pandas-dev/pandas/tree/master/doc/data/air_quality_pm25_long.csv" class="btn btn-dark btn-sm">To raw data</a>
            </div>
        </div>

.. ipython:: python

    air_quality_pm25 = pd.read_csv("data/air_quality_pm25_long.csv",
                                   parse_dates=True)
    air_quality_pm25 = air_quality_pm25[["date.utc", "location",
                                         "parameter", "value"]]
    air_quality_pm25.head()

.. raw:: html

        </li>
    </ul>
    </div>


How to combine data from multiple tables?
-----------------------------------------

Concatenating objects
~~~~~~~~~~~~~~~~~~~~~

.. image:: ../../_static/schemas/08_concat_row.svg
   :align: center

.. raw:: html

    <ul class="task-bullet">
        <li>

I want to combine the measurements of :math:`NO_2` and :math:`PM_{25}`, two tables with a similar structure, in a single table

.. ipython:: python

    air_quality = pd.concat([air_quality_pm25, air_quality_no2], axis=0)
    air_quality.head()

The :func:`~pandas.concat` function performs concatenation operations of multiple
tables along one of the axis (row-wise or column-wise).

.. raw:: html

        </li>
    </ul>

By default concatenation is along axis 0, so the resulting table combines the rows
of the input tables. Let’s check the shape of the original and the
concatenated tables to verify the operation:

.. ipython:: python

    print('Shape of the `air_quality_pm25` table: ', air_quality_pm25.shape)
    print('Shape of the `air_quality_no2` table: ', air_quality_no2.shape)
    print('Shape of the resulting `air_quality` table: ', air_quality.shape)

Hence, the resulting table has 3178 = 1110 + 2068 rows.

.. note::
    The **axis** argument will return in a number of pandas
    methods that can be applied **along an axis**. A ``DataFrame`` has two
    corresponding axes: the first running vertically downwards across rows
    (axis 0), and the second running horizontally across columns (axis 1).
    Most operations like concatenation or summary statistics are by default
    across rows (axis 0), but can be applied across columns as well.

Sorting the table on the datetime information illustrates also the
combination of both tables, with the ``parameter`` column defining the
origin of the table (either ``no2`` from table ``air_quality_no2`` or
``pm25`` from table ``air_quality_pm25``):

.. ipython:: python

    air_quality = air_quality.sort_values("date.utc")
    air_quality.head()

In this specific example, the ``parameter`` column provided by the data
ensures that each of the original tables can be identified. This is not
always the case. the ``concat`` function provides a convenient solution
with the ``keys`` argument, adding an additional (hierarchical) row
index. For example:

.. ipython:: python

    air_quality_ = pd.concat([air_quality_pm25, air_quality_no2],
                             keys=["PM25", "NO2"])

.. ipython:: python

    air_quality_.head()

.. note::
    The existence of multiple row/column indices at the same time
    has not been mentioned within these tutorials. *Hierarchical indexing*
    or *MultiIndex* is an advanced and powerfull pandas feature to analyze
    higher dimensional data.

    Multi-indexing is out of scope for this pandas introduction. For the
    moment, remember that the function ``reset_index`` can be used to
    convert any level of an index to a column, e.g.
    ``air_quality.reset_index(level=0)``

    .. raw:: html

        <div class="d-flex flex-row  gs-torefguide">
            <span class="badge badge-info">To user guide</span>

    Feel free to dive into the world of multi-indexing at the user guide section on :ref:`advanced indexing <advanced>`.

    .. raw:: html

        </div>

.. raw:: html

    <div class="d-flex flex-row gs-torefguide">
        <span class="badge badge-info">To user guide</span>

More options on table concatenation (row and column
wise) and how ``concat`` can be used to define the logic (union or
intersection) of the indexes on the other axes is provided at the section on
:ref:`object concatenation <merging.concat>`.

.. raw:: html

    </div>

Join tables using a common identifier
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. image:: ../../_static/schemas/08_merge_left.svg
   :align: center

.. raw:: html

    <ul class="task-bullet">
        <li>

Add the station coordinates, provided by the stations metadata table, to the corresponding rows in the measurements table.

.. warning::
    The air quality measurement station coordinates are stored in a data
    file ``air_quality_stations.csv``, downloaded using the
    `py-openaq <http://dhhagan.github.io/py-openaq/index.html>`__ package.

.. ipython:: python

    stations_coord = pd.read_csv("data/air_quality_stations.csv")
    stations_coord.head()

.. note::
    The stations used in this example (FR04014, BETR801 and London
    Westminster) are just three entries enlisted in the metadata table. We
    only want to add the coordinates of these three to the measurements
    table, each on the corresponding rows of the ``air_quality`` table.

.. ipython:: python

    air_quality.head()

.. ipython:: python

    air_quality = pd.merge(air_quality, stations_coord,
                           how='left', on='location')
    air_quality.head()

Using the :meth:`~pandas.merge` function, for each of the rows in the
``air_quality`` table, the corresponding coordinates are added from the
``air_quality_stations_coord`` table. Both tables have the column
``location`` in common which is used as a key to combine the
information. By choosing the ``left`` join, only the locations available
in the ``air_quality`` (left) table, i.e. FR04014, BETR801 and London
Westminster, end up in the resulting table. The ``merge`` function
supports multiple join options similar to database-style operations.

.. raw:: html

        </li>
    </ul>

.. raw:: html

    <ul class="task-bullet">
        <li>

Add the parameter full description and name, provided by the parameters metadata table, to the measurements table

.. warning::
    The air quality parameters metadata are stored in a data file
    ``air_quality_parameters.csv``, downloaded using the
    `py-openaq <http://dhhagan.github.io/py-openaq/index.html>`__ package.

.. ipython:: python

    air_quality_parameters = pd.read_csv("data/air_quality_parameters.csv")
    air_quality_parameters.head()

.. ipython:: python

    air_quality = pd.merge(air_quality, air_quality_parameters,
                           how='left', left_on='parameter', right_on='id')
    air_quality.head()

Compared to the previous example, there is no common column name.
However, the ``parameter`` column in the ``air_quality`` table and the
``id`` column in the ``air_quality_parameters_name`` both provide the
measured variable in a common format. The ``left_on`` and ``right_on``
arguments are used here (instead of just ``on``) to make the link
between the two tables.

.. raw:: html

        </li>
    </ul>

.. raw:: html

    <div class="d-flex flex-row gs-torefguide">
        <span class="badge badge-info">To user guide</span>

pandas supports also inner, outer, and right joins.
More information on join/merge of tables is provided in the user guide section on
:ref:`database style merging of tables <merging.join>`. Or have a look at the
:ref:`comparison with SQL<compare_with_sql.join>` page.

.. raw:: html

   </div>

.. raw:: html

    <div class="shadow gs-callout gs-callout-remember">
        <h4>REMEMBER</h4>

-  Multiple tables can be concatenated both column as row wise using
   the ``concat`` function.
-  For database-like merging/joining of tables, use the ``merge``
   function.

.. raw:: html

   </div>

.. raw:: html

    <div class="d-flex flex-row gs-torefguide">
        <span class="badge badge-info">To user guide</span>

See the user guide for a full description of the various :ref:`facilities to combine data tables <merging>`.

.. raw:: html

   </div>
