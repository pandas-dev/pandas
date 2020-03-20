{{ header }}

.. _getting_started:

===============
Getting started
===============

Installation
------------

.. raw:: html

    <div class="container">
        <div class="row">
            <div class="col-lg-6 col-md-6 col-sm-12 col-xs-12 d-flex install-block">
                <div class="card install-card shadow w-100">
                <div class="card-header">
                    Working with conda?
                </div>
                <div class="card-body">
                    <p class="card-text">

pandas is part of the `Anaconda <http://docs.continuum.io/anaconda/>`__ distribution and can be
installed with Anaconda or Miniconda:

.. raw:: html

                    </p>
                </div>
                <div class="card-footer text-muted">

.. code-block:: bash

   conda install pandas

.. raw:: html

                </div>
                </div>
            </div>
            <div class="col-lg-6 col-md-6 col-sm-12 col-xs-12 d-flex install-block">
                <div class="card install-card shadow w-100">
                <div class="card-header">
                    Prefer pip?
                </div>
                <div class="card-body">
                    <p class="card-text">

pandas can be installed via pip from `PyPI <https://pypi.org/project/pandas>`__.

.. raw:: html

                    </p>
                </div>
                <div class="card-footer text-muted">

.. code-block:: bash

   pip install pandas

.. raw:: html

                </div>
                </div>
            </div>
            <div class="col-12 d-flex install-block">
                <div class="card install-card shadow w-100">
                <div class="card-header">
                    In-depth instructions?
                </div>
                <div class="card-body">
                    <p class="card-text">Installing a specific version?
                      Installing from source?
                      Check the advanced installation page.</p>

.. container:: custom-button

    :ref:`Learn more <install>`

.. raw:: html

                </div>
                </div>
            </div>
        </div>
    </div>

.. _gentle_intro:

Intro to pandas
---------------

.. raw:: html

    <div class="container">
    <div id="accordion" class="shadow tutorial-accordion">

        <div class="card tutorial-card">
            <div class="card-header collapsed card-link" data-toggle="collapse" data-target="#collapseOne">
                <div class="d-flex flex-row tutorial-card-header-1">
                    <div class="d-flex flex-row tutorial-card-header-2">
                        <button class="btn btn-dark btn-sm"></button>
                        What kind of data does pandas handle?
                    </div>
                    <span class="badge gs-badge-link">

:ref:`Straight to tutorial...<10min_tut_01_tableoriented>`

.. raw:: html

                    </span>
                </div>
            </div>
            <div id="collapseOne" class="collapse" data-parent="#accordion">
                <div class="card-body">

When working with tabular data, such as data stored in spreadsheets or databases, pandas is the right tool for you. pandas will help you
to explore, clean and process your data. In pandas, a data table is called a :class:`DataFrame`.

.. image:: ../_static/schemas/01_table_dataframe.svg
   :align: center

.. raw:: html

                    <div class="d-flex flex-row">
                        <span class="badge gs-badge-link">

:ref:`To introduction tutorial <10min_tut_01_tableoriented>`

.. raw:: html

                        </span>
                        <span class="badge gs-badge-link">

:ref:`To user guide <dsintro>`

.. raw:: html

                        </span>
                    </div>
                </div>
            </div>
        </div>

        <div class="card tutorial-card">
            <div class="card-header collapsed card-link" data-toggle="collapse" data-target="#collapseTwo">
                <div class="d-flex flex-row tutorial-card-header-1">
                    <div class="d-flex flex-row tutorial-card-header-2">
                        <button class="btn btn-dark btn-sm"></button>
                        How do I read and write tabular data?
                    </div>
                    <span class="badge gs-badge-link">

:ref:`Straight to tutorial...<10min_tut_02_read_write>`

.. raw:: html

                    </span>
                </div>
            </div>
            <div id="collapseTwo" class="collapse" data-parent="#accordion">
                <div class="card-body">

pandas supports the integration with many file formats or data sources out of the box (csv, excel, sql, json, parquet,…). Importing data from each of these
data sources is provided by function with the prefix ``read_*``. Similarly, the ``to_*`` methods are used to store data.

.. image:: ../_static/schemas/02_io_readwrite.svg
   :align: center

.. raw:: html

                    <div class="d-flex flex-row">
                        <span class="badge gs-badge-link">

:ref:`To introduction tutorial <10min_tut_02_read_write>`

.. raw:: html

                        </span>
                        <span class="badge gs-badge-link">

:ref:`To user guide <io>`

.. raw:: html

                        </span>
                    </div>
                </div>
            </div>
        </div>

        <div class="card tutorial-card">
            <div class="card-header collapsed card-link" data-toggle="collapse" data-target="#collapseThree">
                <div class="d-flex flex-row tutorial-card-header-1">
                    <div class="d-flex flex-row tutorial-card-header-2">
                        <button class="btn btn-dark btn-sm"></button>
                        How do I select a subset of a table?
                    </div>
                    <span class="badge gs-badge-link">

:ref:`Straight to tutorial...<10min_tut_03_subset>`

.. raw:: html

                    </span>
                </div>
            </div>
            <div id="collapseThree" class="collapse" data-parent="#accordion">
                <div class="card-body">

Selecting or filtering specific rows and/or columns? Filtering the data on a condition? Methods for slicing, selecting, and extracting the
data you need are available in pandas.

.. image:: ../_static/schemas/03_subset_columns_rows.svg
   :align: center

.. raw:: html

                    <div class="d-flex flex-row">
                        <span class="badge gs-badge-link">

:ref:`To introduction tutorial <10min_tut_03_subset>`

.. raw:: html

                        </span>
                        <span class="badge gs-badge-link">

:ref:`To user guide <indexing>`

.. raw:: html

                        </span>
                    </div>
                </div>
            </div>
        </div>

        <div class="card tutorial-card">
            <div class="card-header collapsed card-link" data-toggle="collapse" data-target="#collapseFour">
                <div class="d-flex flex-row tutorial-card-header-1">
                    <div class="d-flex flex-row tutorial-card-header-2">
                        <button class="btn btn-dark btn-sm"></button>
                        How to create plots in pandas?
                    </div>
                    <span class="badge gs-badge-link">

:ref:`Straight to tutorial...<10min_tut_04_plotting>`

.. raw:: html

                    </span>
                </div>
            </div>
            <div id="collapseFour" class="collapse" data-parent="#accordion">
                <div class="card-body">

pandas provides plotting your data out of the box, using the power of Matplotlib. You can pick the plot type (scatter, bar, boxplot,...)
corresponding to your data.

.. image:: ../_static/schemas/04_plot_overview.svg
   :align: center

.. raw:: html

                    <div class="d-flex flex-row">
                        <span class="badge gs-badge-link">

:ref:`To introduction tutorial <10min_tut_04_plotting>`

.. raw:: html

                        </span>
                        <span class="badge gs-badge-link">

:ref:`To user guide <visualization>`

.. raw:: html

                        </span>
                    </div>
                </div>
            </div>
        </div>

        <div class="card tutorial-card">
            <div class="card-header collapsed card-link" data-toggle="collapse" data-target="#collapseFive">
                <div class="d-flex flex-row tutorial-card-header-1">
                    <div class="d-flex flex-row tutorial-card-header-2">
                        <button class="btn btn-dark btn-sm"></button>
                        How to create new columns derived from existing columns?
                    </div>
                    <span class="badge gs-badge-link">

:ref:`Straight to tutorial...<10min_tut_05_columns>`

.. raw:: html

                    </span>
                </div>
            </div>
            <div id="collapseFive" class="collapse" data-parent="#accordion">
                <div class="card-body">

There is no need to loop over all rows of your data table to do calculations. Data manipulations on a column work elementwise.
Adding a column to a :class:`DataFrame` based on existing data in other columns is straightforward.

.. image:: ../_static/schemas/05_newcolumn_2.svg
   :align: center

.. raw:: html

                    <div class="d-flex flex-row">
                        <span class="badge gs-badge-link">

:ref:`To introduction tutorial <10min_tut_05_columns>`

.. raw:: html

                        </span>
                        <span class="badge gs-badge-link">

:ref:`To user guide <basics.dataframe.sel_add_del>`

.. raw:: html

                        </span>
                    </div>
                </div>
            </div>
        </div>

        <div class="card tutorial-card">
            <div class="card-header collapsed card-link" data-toggle="collapse" data-target="#collapseSix">
                <div class="d-flex flex-row tutorial-card-header-1">
                    <div class="d-flex flex-row tutorial-card-header-2">
                        <button class="btn btn-dark btn-sm"></button>
                        How to calculate summary statistics?
                    </div>
                    <span class="badge gs-badge-link">

:ref:`Straight to tutorial...<10min_tut_06_stats>`

.. raw:: html

                    </span>
                </div>
            </div>
            <div id="collapseSix" class="collapse" data-parent="#accordion">
                <div class="card-body">

Basic statistics (mean, median, min, max, counts...) are easily calculable. These or custom aggregations can be applied on the entire
data set, a sliding window of the data or grouped by categories. The latter is also known as the split-apply-combine approach.

.. image:: ../_static/schemas/06_groupby.svg
   :align: center

.. raw:: html

                    <div class="d-flex flex-row">
                        <span class="badge gs-badge-link">

:ref:`To introduction tutorial <10min_tut_06_stats>`

.. raw:: html

                        </span>
                        <span class="badge gs-badge-link">

:ref:`To user guide <groupby>`

.. raw:: html

                        </span>
                    </div>
                </div>
            </div>
        </div>

        <div class="card tutorial-card">
            <div class="card-header collapsed card-link" data-toggle="collapse" data-target="#collapseSeven">
                <div class="d-flex flex-row tutorial-card-header-1">
                    <div class="d-flex flex-row tutorial-card-header-2">
                        <button class="btn btn-dark btn-sm"></button>
                        How to reshape the layout of tables?
                    </div>
                    <span class="badge gs-badge-link">

:ref:`Straight to tutorial...<10min_tut_07_reshape>`

.. raw:: html

                    </span>
                </div>
            </div>
            <div id="collapseSeven" class="collapse" data-parent="#accordion">
                <div class="card-body">

Change the structure of your data table in multiple ways. You can :func:`~pandas.melt` your data table from wide to long/tidy form or :func:`~pandas.pivot`
from long to wide format. With aggregations built-in, a pivot table is created with a sinlge command.

.. image:: ../_static/schemas/07_melt.svg
   :align: center

.. raw:: html

                    <div class="d-flex flex-row">
                        <span class="badge gs-badge-link">

:ref:`To introduction tutorial <10min_tut_07_reshape>`

.. raw:: html

                        </span>
                        <span class="badge gs-badge-link">

:ref:`To user guide <reshaping>`

.. raw:: html

                        </span>
                    </div>
                </div>
            </div>
        </div>

        <div class="card tutorial-card">
            <div class="card-header collapsed card-link" data-toggle="collapse" data-target="#collapseEight">
                <div class="d-flex flex-row tutorial-card-header-1">
                    <div class="d-flex flex-row tutorial-card-header-2">
                        <button class="btn btn-dark btn-sm"></button>
                        How to combine data from multiple tables?
                    </div>
                    <span class="badge gs-badge-link">

:ref:`Straight to tutorial...<10min_tut_08_combine>`

.. raw:: html

                    </span>
                </div>
            </div>
            <div id="collapseEight" class="collapse" data-parent="#accordion">
                <div class="card-body">

Multiple tables can be concatenated both column wise as row wise and database-like join/merge operations are provided to combine multiple tables of data.

.. image:: ../_static/schemas/08_concat_row.svg
   :align: center

.. raw:: html

                    <div class="d-flex flex-row">
                        <span class="badge gs-badge-link">

:ref:`To introduction tutorial <10min_tut_08_combine>`

.. raw:: html

                        </span>
                        <span class="badge gs-badge-link">

:ref:`To user guide <merging>`

.. raw:: html

                        </span>
                    </div>
                </div>
            </div>
        </div>

        <div class="card tutorial-card">
            <div class="card-header collapsed card-link" data-toggle="collapse" data-target="#collapseNine">
                <div class="d-flex flex-row tutorial-card-header-1">
                    <div class="d-flex flex-row tutorial-card-header-2">
                        <button class="btn btn-dark btn-sm"></button>
                        How to handle time series data?
                    </div>
                    <span class="badge gs-badge-link">

:ref:`Straight to tutorial...<10min_tut_09_timeseries>`

.. raw:: html

                    </span>
                </div>
            </div>
            <div id="collapseNine" class="collapse" data-parent="#accordion">
                <div class="card-body">

pandas has great support for time series and has an extensive set of tools for working with dates, times, and time-indexed data.

.. raw:: html

                    <div class="d-flex flex-row">
                        <span class="badge gs-badge-link">

:ref:`To introduction tutorial <10min_tut_09_timeseries>`

.. raw:: html

                        </span>
                        <span class="badge gs-badge-link">

:ref:`To user guide <timeseries>`

.. raw:: html

                        </span>
                    </div>
                </div>
            </div>
        </div>

        <div class="card tutorial-card">
            <div class="card-header collapsed card-link" data-toggle="collapse" data-target="#collapseTen">
                <div class="d-flex flex-row tutorial-card-header-1">
                    <div class="d-flex flex-row tutorial-card-header-2">
                        <button class="btn btn-dark btn-sm"></button>
                        How to manipulate textual data?
                    </div>
                    <span class="badge gs-badge-link">

:ref:`Straight to tutorial...<10min_tut_10_text>`

.. raw:: html

                    </span>
                </div>
            </div>
            <div id="collapseTen" class="collapse" data-parent="#accordion">
                <div class="card-body">

Data sets do not only contain numerical data. pandas provides a wide range of functions to cleaning textual data and extract useful information from it.

.. raw:: html

                    <div class="d-flex flex-row">
                        <span class="badge gs-badge-link">

:ref:`To introduction tutorial <10min_tut_10_text>`

.. raw:: html

                        </span>
                        <span class="badge gs-badge-link">

:ref:`To user guide <timeseries>`

.. raw:: html

                        </span>
                    </div>
                </div>
            </div>
        </div>

    </div>
    </div>


.. _comingfrom:

Coming from...
--------------

Are you familiar with other software for manipulating tablular data? Learn
the pandas-equivalent operations compared to software you already know:

.. raw:: html

    <div class="container">
        <div class="row">
            <div class="col-lg-6 col-md-6 col-sm-6 col-xs-12 d-flex">
                <div class="card text-center intro-card shadow">
                <img src="../_static/logo_r.svg" class="card-img-top" alt="R project logo" height="72">
                <div class="card-body flex-fill">
                    <p class="card-text">The <a href="https://www.r-project.org/">R programming language</a> provides the <code>dataframe</code> data structure and multiple packages,
                        such as <a href="https://www.tidyverse.org/">tidyverse</a> use and extend <code>data.frame</code>s for convenient data handling
                        functionalities similar to pandas.</p>

.. container:: custom-button

    :ref:`Learn more <compare_with_r>`

.. raw:: html

                </div>
                </div>
            </div>
            <div class="col-lg-6 col-md-6 col-sm-6 col-xs-12 d-flex">
                <div class="card text-center intro-card shadow">
                <img src="../_static/logo_sql.svg" class="card-img-top" alt="SQL logo" height="72">
                <div class="card-body flex-fill">
                    <p class="card-text">Already familiar to <code>SELECT</code>, <code>GROUP BY</code>, <code>JOIN</code>, etc.?
                    Most of these SQL manipulations do have equivalents in pandas.</p>

.. container:: custom-button

    :ref:`Learn more <compare_with_sql>`

.. raw:: html

                    </div>
                    </div>
                </div>
                <div class="col-lg-6 col-md-6 col-sm-6 col-xs-12 d-flex">
                    <div class="card text-center intro-card shadow">
                    <img src="../_static/logo_stata.svg" class="card-img-top" alt="STATA logo" height="52">
                    <div class="card-body flex-fill">
                        <p class="card-text">The <code>data set</code> included in the
                            <a href="https://en.wikipedia.org/wiki/Stata">STATA</a> statistical software suite corresponds
                            to the pandas <code>dataframe</code>. Many of the operations known from STATA have an equivalent
                            in pandas.</p>

.. container:: custom-button

    :ref:`Learn more <compare_with_stata>`

.. raw:: html

                    </div>
                    </div>
                </div>
                <div class="col-lg-6 col-md-6 col-sm-6 col-xs-12 d-flex">
                    <div class="card text-center intro-card shadow">
                    <img src="../_static/logo_sas.svg" class="card-img-top" alt="SAS logo" height="52">
                    <div class="card-body flex-fill">
                        <p class="card-text">The  <a href="https://en.wikipedia.org/wiki/SAS_(software)">SAS</a> statistical software suite
                            also provides the <code>data set</code> corresponding to the pandas <code>dataframe</code>.
                            Also SAS vectorized operations, filtering, string processing operations, and more have similar
                            functions in pandas.</p>

.. container:: custom-button

    :ref:`Learn more <compare_with_sas>`

.. raw:: html

                    </div>
                    </div>
                </div>
        </div>
    </div>

Tutorials
---------

For a quick overview of pandas functionality, see :ref:`10 Minutes to pandas<10min>`.

You can also reference the pandas `cheat sheet <https://pandas.pydata.org/Pandas_Cheat_Sheet.pdf>`_
for a succinct guide for manipulating data with pandas.

The community produces a wide variety of tutorials available online. Some of the
material is enlisted in the community contributed :ref:`communitytutorials`.


.. If you update this toctree, also update the manual toctree in the
   main index.rst.template

.. toctree::
    :maxdepth: 2
    :hidden:

    install
    overview
    intro_tutorials/index
    comparison/index
    tutorials
