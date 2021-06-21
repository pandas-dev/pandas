.. _10min_tut_04_plotting:

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

.. include:: includes/air_quality_no2.rst

.. ipython:: python

    air_quality = pd.read_csv("data/air_quality_no2.csv", index_col=0, parse_dates=True)
    air_quality.head()

.. note::
    The usage of the ``index_col`` and ``parse_dates`` parameters of the ``read_csv`` function to define the first (0th) column as
    index of the resulting ``DataFrame`` and convert the dates in the column to :class:`Timestamp` objects, respectively.

.. raw:: html

        </li>
    </ul>
    </div>

How to create plots in pandas?
------------------------------

.. image:: ../../_static/schemas/04_plot_overview.svg
   :align: center

.. raw:: html

    <ul class="task-bullet">
        <li>

I want a quick visual check of the data.

.. ipython:: python

    @savefig 04_airqual_quick.png
    air_quality.plot()

With a ``DataFrame``, pandas creates by default one line plot for each of
the columns with numeric data.

.. raw:: html

        </li>
    </ul>

.. raw:: html

    <ul class="task-bullet">
        <li>

I want to plot only the columns of the data table with the data from Paris.

.. ipython:: python

    @savefig 04_airqual_paris.png
    air_quality["station_paris"].plot()

To plot a specific column, use the selection method of the
:ref:`subset data tutorial <10min_tut_03_subset>` in combination with the :meth:`~DataFrame.plot`
method. Hence, the :meth:`~DataFrame.plot` method works on both ``Series`` and
``DataFrame``.

.. raw:: html

        </li>
    </ul>

.. raw:: html

    <ul class="task-bullet">
        <li>

I want to visually compare the :math:`N0_2` values measured in London versus Paris.

.. ipython:: python

    @savefig 04_airqual_scatter.png
    air_quality.plot.scatter(x="station_london", y="station_paris", alpha=0.5)

.. raw:: html

        </li>
    </ul>

Apart from the default ``line`` plot when using the ``plot`` function, a
number of alternatives are available to plot data. Let’s use some
standard Python to get an overview of the available plot methods:

.. ipython:: python

    [
        method_name
        for method_name in dir(air_quality.plot)
        if not method_name.startswith("_")
    ]

.. note::
    In many development environments as well as IPython and
    Jupyter Notebook, use the TAB button to get an overview of the available
    methods, for example ``air_quality.plot.`` + TAB.

One of the options is :meth:`DataFrame.plot.box`, which refers to a
`boxplot <https://en.wikipedia.org/wiki/Box_plot>`__. The ``box``
method is applicable on the air quality example data:

.. ipython:: python

    @savefig 04_airqual_boxplot.png
    air_quality.plot.box()

.. raw:: html

    <div class="d-flex flex-row gs-torefguide">
        <span class="badge badge-info">To user guide</span>

For an introduction to plots other than the default line plot, see the user guide section about :ref:`supported plot styles <visualization.other>`.

.. raw:: html

   </div>

.. raw:: html

    <ul class="task-bullet">
        <li>

I want each of the columns in a separate subplot.

.. ipython:: python

    @savefig 04_airqual_area_subplot.png
    axs = air_quality.plot.area(figsize=(12, 4), subplots=True)

Separate subplots for each of the data columns are supported by the ``subplots`` argument
of the ``plot`` functions. The builtin options available in each of the pandas plot
functions are worth reviewing.

.. raw:: html

        </li>
    </ul>

.. raw:: html

    <div class="d-flex flex-row gs-torefguide">
        <span class="badge badge-info">To user guide</span>

Some more formatting options are explained in the user guide section on :ref:`plot formatting <visualization.formatting>`.

.. raw:: html

   </div>

.. raw:: html

    <ul class="task-bullet">
        <li>

I want to further customize, extend or save the resulting plot.

.. ipython:: python

    fig, axs = plt.subplots(figsize=(12, 4))
    air_quality.plot.area(ax=axs)
    @savefig 04_airqual_customized.png
    axs.set_ylabel("NO$_2$ concentration")
    fig.savefig("no2_concentrations.png")

.. ipython:: python
   :suppress:

   import os

   os.remove("no2_concentrations.png")

.. raw:: html

        </li>
    </ul>

Each of the plot objects created by pandas is a
`matplotlib <https://matplotlib.org/>`__ object. As Matplotlib provides
plenty of options to customize plots, making the link between pandas and
Matplotlib explicit enables all the power of matplotlib to the plot.
This strategy is applied in the previous example:

::

   fig, axs = plt.subplots(figsize=(12, 4))        # Create an empty matplotlib Figure and Axes
   air_quality.plot.area(ax=axs)                   # Use pandas to put the area plot on the prepared Figure/Axes
   axs.set_ylabel("NO$_2$ concentration")          # Do any matplotlib customization you like
   fig.savefig("no2_concentrations.png")           # Save the Figure/Axes using the existing matplotlib method.

.. raw:: html

    <div class="shadow gs-callout gs-callout-remember">
        <h4>REMEMBER</h4>

-  The ``.plot.*`` methods are applicable on both Series and DataFrames
-  By default, each of the columns is plotted as a different element
   (line, boxplot,…)
-  Any plot created by pandas is a Matplotlib object.

.. raw:: html

   </div>

.. raw:: html

    <div class="d-flex flex-row gs-torefguide">
        <span class="badge badge-info">To user guide</span>

A full overview of plotting in pandas is provided in the :ref:`visualization pages <visualization>`.

.. raw:: html

   </div>
