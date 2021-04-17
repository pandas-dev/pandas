{{ header }}

.. _getting_started:

===============
Getting started
===============

Installation
------------

.. panels::
    :card: + install-card
    :column: col-lg-6 col-md-6 col-sm-12 col-xs-12 p-3

    Working with conda?
    ^^^^^^^^^^^^^^^^^^^

    pandas is part of the `Anaconda <https://docs.continuum.io/anaconda/>`__
    distribution and can be installed with Anaconda or Miniconda:

    ++++++++++++++++++++++

    .. code-block:: bash

        conda install pandas

    ---

    Prefer pip?
    ^^^^^^^^^^^

    pandas can be installed via pip from `PyPI <https://pypi.org/project/pandas>`__.

    ++++

    .. code-block:: bash

        pip install pandas

    ---
    :column: col-12 p-3

    In-depth instructions?
    ^^^^^^^^^^^^^^^^^^^^^^

    Installing a specific version? Installing from source? Check the advanced
    installation page.

    .. link-button:: ./install.html
        :type: url
        :text: Learn more
        :classes: btn-secondary stretched-link


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
to explore, clean, and process your data. In pandas, a data table is called a :class:`DataFrame`.

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
data set, a sliding window of the data, or grouped by categories. The latter is also known as the split-apply-combine approach.

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
from long to wide format. With aggregations built-in, a pivot table is created with a single command.

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

Multiple tables can be concatenated both column wise and row wise as database-like join/merge operations are provided to combine multiple tables of data.

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

Data sets do not only contain numerical data. pandas provides a wide range of functions to clean textual data and extract useful information from it.

.. raw:: html

                    <div class="d-flex flex-row">
                        <span class="badge gs-badge-link">

:ref:`To introduction tutorial <10min_tut_10_text>`

.. raw:: html

                        </span>
                        <span class="badge gs-badge-link">

:ref:`To user guide <text>`

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

.. panels::
    :card: + comparison-card text-center shadow
    :column: col-lg-6 col-md-6 col-sm-6 col-xs-12 d-flex

    ---
    :card: + comparison-card-r
    :img-top: ../_static/logo_r.svg

    The `R programming language <https://www.r-project.org/>`__ provides the
    ``data.frame`` data structure and multiple packages, such as
    `tidyverse <https://www.tidyverse.org>`__ use and extend ``data.frame``
    for convenient data handling functionalities similar to pandas.

    +++

    .. link-button:: compare_with_r
        :type: ref
        :text: Learn more
        :classes: btn-secondary stretched-link


    ---
    :card: + comparison-card-sql
    :img-top: ../_static/logo_sql.svg

    Already familiar to ``SELECT``, ``GROUP BY``, ``JOIN``, etc.?
    Most of these SQL manipulations do have equivalents in pandas.

    +++

    .. link-button:: compare_with_sql
        :type: ref
        :text: Learn more
        :classes: btn-secondary stretched-link


    ---
    :card: + comparison-card-stata
    :img-top: ../_static/logo_stata.svg

    The ``data set`` included in the `STATA <https://en.wikipedia.org/wiki/Stata>`__
    statistical software suite corresponds to the pandas ``DataFrame``.
    Many of the operations known from STATA have an equivalent in pandas.

    +++

    .. link-button:: compare_with_stata
        :type: ref
        :text: Learn more
        :classes: btn-secondary stretched-link


    ---
    :card: + comparison-card-excel
    :img-top: ../_static/spreadsheets/logo_excel.svg

    Users of `Excel <https://en.wikipedia.org/wiki/Microsoft_Excel>`__
    or other spreadsheet programs will find that many of the concepts are
    transferrable to pandas.

    +++

    .. link-button:: compare_with_spreadsheets
        :type: ref
        :text: Learn more
        :classes: btn-secondary stretched-link


    ---
    :card: + comparison-card-sas
    :img-top: ../_static/logo_sas.svg

    The `SAS <https://en.wikipedia.org/wiki/SAS_(software)>`__ statistical software suite
    also provides the ``data set`` corresponding to the pandas ``DataFrame``.
    Also SAS vectorized operations, filtering, string processing operations,
    and more have similar functions in pandas.

    +++

    .. link-button:: compare_with_sas
        :type: ref
        :text: Learn more
        :classes: btn-secondary stretched-link


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
