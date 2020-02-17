.. _10min_tut_09_timeseries:

{{ header }}

.. ipython:: python

    import pandas as pd
    import matplotlib.pyplot as plt

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
                    <span class="badge badge-dark">Air quality data</span>
                </div>
                <div class="collapse" id="collapsedata">
                    <div class="card-body">
                        <p class="card-text">

For this tutorial, air quality data about :math:`NO_2` and Particulate
matter less than 2.5 micrometers is used, made available by
`openaq <https://openaq.org>`__ and downloaded using the
`py-openaq <http://dhhagan.github.io/py-openaq/index.html>`__ package.
The ``air_quality_no2_long.csv"`` data set provides :math:`NO_2` values
for the measurement stations *FR04014*, *BETR801* and *London
Westminster* in respectively Paris, Antwerp and London.

.. raw:: html

                        </p>
                    <a href="https://github.com/pandas-dev/pandas/tree/master/doc/data/air_quality_no2_long.csv" class="btn btn-dark btn-sm">To raw data</a>
                </div>
            </div>

.. ipython:: python

    air_quality = pd.read_csv("data/air_quality_no2_long.csv")
    air_quality = air_quality.rename(columns={"date.utc": "datetime"})
    air_quality.head()

.. ipython:: python

    air_quality.city.unique()

.. raw:: html

        </li>
    </ul>
    </div>

How to handle time series data with ease?
-----------------------------------------

Using pandas datetime properties
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. raw:: html

    <ul class="task-bullet">
        <li>

I want to work with the dates in the column ``datetime`` as datetime objects instead of plain text

.. ipython:: python

    air_quality["datetime"] = pd.to_datetime(air_quality["datetime"])
    air_quality["datetime"]

Initially, the values in ``datetime`` are character strings and do not
provide any datetime operations (e.g. extract the year, day of the
week,…). By applying the ``to_datetime`` function, pandas interprets the
strings and convert these to datetime (i.e. ``datetime64[ns, UTC]``)
objects. In pandas we call these datetime objects similar to
``datetime.datetime`` from the standard library a :class:`pandas.Timestamp`.

.. raw:: html

        </li>
    </ul>

.. note::
    As many data sets do contain datetime information in one of
    the columns, pandas input function like :func:`pandas.read_csv` and :func:`pandas.read_json`
    can do the transformation to dates when reading the data using the
    ``parse_dates`` parameter with a list of the columns to read as
    Timestamp:

    ::

        pd.read_csv("../data/air_quality_no2_long.csv", parse_dates=["datetime"])

Why are these :class:`pandas.Timestamp` objects useful. Let’s illustrate the added
value with some example cases.

   What is the start and end date of the time series data set working
   with?

.. ipython:: python

    air_quality["datetime"].min(), air_quality["datetime"].max()

Using :class:`pandas.Timestamp` for datetimes enable us to calculate with date
information and make them comparable. Hence, we can use this to get the
length of our time series:

.. ipython:: python

    air_quality["datetime"].max() - air_quality["datetime"].min()

The result is a :class:`pandas.Timedelta` object, similar to ``datetime.timedelta``
from the standard Python library and defining a time duration.

.. raw:: html

    <div class="d-flex flex-row gs-torefguide">
        <span class="badge badge-info">To user guide</span>

The different time concepts supported by pandas are explained in the user guide section on :ref:`time related concepts <timeseries.overview>`.

.. raw:: html

    </div>

.. raw:: html

    <ul class="task-bullet">
        <li>

I want to add a new column to the ``DataFrame`` containing only the month of the measurement

.. ipython:: python

    air_quality["month"] = air_quality["datetime"].dt.month
    air_quality.head()

By using ``Timestamp`` objects for dates, a lot of time-related
properties are provided by pandas. For example the ``month``, but also
``year``, ``weekofyear``, ``quarter``,… All of these properties are
accessible by the ``dt`` accessor.

.. raw:: html

        </li>
    </ul>

.. raw:: html

    <div class="d-flex flex-row gs-torefguide">
        <span class="badge badge-info">To user guide</span>

An overview of the existing date properties is given in the
:ref:`time and date components overview table <timeseries.components>`. More details about the ``dt`` accessor
to return datetime like properties is explained in a dedicated section on the  :ref:`dt accessor <basics.dt_accessors>`.

.. raw:: html

    </div>

.. raw:: html

    <ul class="task-bullet">
        <li>

What is the average :math:`NO_2` concentration for each day of the week for each of the measurement locations?

.. ipython:: python

    air_quality.groupby(
        [air_quality["datetime"].dt.weekday, "location"])["value"].mean()

Remember the split-apply-combine pattern provided by ``groupby`` from the
:ref:`tutorial on statistics calculation <10min_tut_06_stats>`?
Here, we want to calculate a given statistic (e.g. mean :math:`NO_2`)
**for each weekday** and **for each measurement location**. To group on
weekdays, we use the datetime property ``weekday`` (with Monday=0 and
Sunday=6) of pandas ``Timestamp``, which is also accessible by the
``dt`` accessor. The grouping on both locations and weekdays can be done
to split the calculation of the mean on each of these combinations.

.. danger::
    As we are working with a very short time series in these
    examples, the analysis does not provide a long-term representative
    result!

.. raw:: html

        </li>
    </ul>

.. raw:: html

    <ul class="task-bullet">
        <li>

Plot the typical :math:`NO_2` pattern during the day of our time series of all stations together. In other words, what is the average value for each hour of the day?

.. ipython:: python

    fig, axs = plt.subplots(figsize=(12, 4))
    air_quality.groupby(
        air_quality["datetime"].dt.hour)["value"].mean().plot(kind='bar',
                                                              rot=0,
                                                              ax=axs)
    plt.xlabel("Hour of the day");  # custom x label using matplotlib
    @savefig 09_bar_chart.png
    plt.ylabel("$NO_2 (µg/m^3)$");

Similar to the previous case, we want to calculate a given statistic
(e.g. mean :math:`NO_2`) **for each hour of the day** and we can use the
split-apply-combine approach again. For this case, the datetime property ``hour``
of pandas ``Timestamp``, which is also accessible by the ``dt`` accessor.

.. raw:: html

        </li>
    </ul>

Datetime as index
~~~~~~~~~~~~~~~~~

In the :ref:`tutorial on reshaping <10min_tut_07_reshape>`,
:meth:`~pandas.pivot` was introduced to reshape the data table with each of the
measurements locations as a separate column:

.. ipython:: python

    no_2 = air_quality.pivot(index="datetime", columns="location", values="value")
    no_2.head()

.. note::
    By pivoting the data, the datetime information became the
    index of the table. In general, setting a column as an index can be
    achieved by the ``set_index`` function.

Working with a datetime index (i.e. ``DatetimeIndex``) provides powerful
functionalities. For example, we do not need the ``dt`` accessor to get
the time series properties, but have these properties available on the
index directly:

.. ipython:: python

    no_2.index.year, no_2.index.weekday

Some other advantages are the convenient subsetting of time period or
the adapted time scale on plots. Let’s apply this on our data.

.. raw:: html

    <ul class="task-bullet">
        <li>

Create a plot of the :math:`NO_2` values in the different stations from the 20th of May till the end of 21st of May

.. ipython:: python
    :okwarning:

    @savefig 09_time_section.png
    no_2["2019-05-20":"2019-05-21"].plot();

By providing a **string that parses to a datetime**, a specific subset of the data can be selected on a ``DatetimeIndex``.

.. raw:: html

        </li>
    </ul>

.. raw:: html

    <div class="d-flex flex-row gs-torefguide">
        <span class="badge badge-info">To user guide</span>

More information on the ``DatetimeIndex`` and the slicing by using strings is provided in the section on :ref:`time series indexing <timeseries.datetimeindex>`.

.. raw:: html

    </div>

Resample a time series to another frequency
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. raw:: html

    <ul class="task-bullet">
        <li>

Aggregate the current hourly time series values to the monthly maximum value in each of the stations.

.. ipython:: python

    monthly_max = no_2.resample("M").max()
    monthly_max

A very powerful method on time series data with a datetime index, is the
ability to :meth:`~Series.resample` time series to another frequency (e.g.,
converting secondly data into 5-minutely data).

.. raw:: html

        </li>
    </ul>

The :meth:`~Series.resample` method is similar to a groupby operation:

-  it provides a time-based grouping, by using a string (e.g. ``M``,
   ``5H``,…) that defines the target frequency
-  it requires an aggregation function such as ``mean``, ``max``,…

.. raw:: html

    <div class="d-flex flex-row gs-torefguide">
        <span class="badge badge-info">To user guide</span>

An overview of the aliases used to define time series frequencies is given in the :ref:`offset aliases overview table <timeseries.offset_aliases>`.

.. raw:: html

    </div>

When defined, the frequency of the time series is provided by the
``freq`` attribute:

.. ipython:: python

    monthly_max.index.freq

.. raw:: html

    <ul class="task-bullet">
        <li>

Make a plot of the daily median :math:`NO_2` value in each of the stations.

.. ipython:: python
    :okwarning:

    @savefig 09_resample_mean.png
    no_2.resample("D").mean().plot(style="-o", figsize=(10, 5));

.. raw:: html

        </li>
    </ul>

.. raw:: html

    <div class="d-flex flex-row gs-torefguide">
        <span class="badge badge-info">To user guide</span>

More details on the power of time series ``resampling`` is provided in the user gudie section on :ref:`resampling <timeseries.resampling>`.

.. raw:: html

    </div>

.. raw:: html

    <div class="shadow gs-callout gs-callout-remember">
        <h4>REMEMBER</h4>

-  Valid date strings can be converted to datetime objects using
   ``to_datetime`` function or as part of read functions.
-  Datetime objects in pandas supports calculations, logical operations
   and convenient date-related properties using the ``dt`` accessor.
-  A ``DatetimeIndex`` contains these date-related properties and
   supports convenient slicing.
-  ``Resample`` is a powerful method to change the frequency of a time
   series.

.. raw:: html

   </div>

.. raw:: html

    <div class="d-flex flex-row gs-torefguide">
        <span class="badge badge-info">To user guide</span>

A full overview on time series is given in the pages on :ref:`time series and date functionality <timeseries>`.

.. raw:: html

   </div>