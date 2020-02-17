{{ header }}

.. _getting_started:

===============
Getting started
===============

Installation
------------

Before you can use pandas, you’ll need to get it installed.

.. container:: container

    .. container:: row

        .. container:: col-lg-6 col-md-6 col-sm-12 col-xs-12 d-flex install-block

            .. container:: card intro-card shadow w-100

                .. container:: card-body

                    .. container:: card-title

                        Working with conda?

                    Pandas is part of the `Anaconda <http://docs.continuum.io/anaconda/>`_ distribution and can be
                    installed with Anaconda or Miniconda:

                .. container:: card-footer text-muted

                    .. code-block:: bash

                        conda install pandas

        .. container:: col-lg-6 col-md-6 col-sm-12 col-xs-12 d-flex install-block

            .. container:: card intro-card shadow w-100

                .. container:: card-body

                    .. container:: card-title

                        Prefer pip?

                    Pandas can be installed via pip from `PyPI <https://pypi.org/project/pandas>`__.

                .. container:: card-footer text-muted

                    .. code-block:: bash

                        pip install pandas

        .. container:: col-12 d-flex install-block

            .. container:: card intro-card shadow w-100

                .. container:: card-body

                    .. container:: card-title

                        In-depth instructions?

                    Installing a specific version? Installing from source? Check the advanced installation page.

                    .. container:: custom-button

                        :ref:`Learn more <install>`

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
                        What kind of data does Pandas handle?
                    </div>
                    <span class="badge gs-badge-link">

:ref:`Straight to tutorial...<10min_tut_01_tableoriented>`

.. raw:: html

                    </span>
                </div>
            </div>
            <div id="collapseOne" class="collapse" data-parent="#accordion">
                <div class="card-body">

When working with tabular data, such as data stored in spreadsheets or databases, Pandas is the right tool for you. Pandas will help you
to explore, clean and process your data. In Pandas, a data table is called a :class:`DataFrame`.

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

Pandas supports the integration with many file formats or data sources out of the box (csv, excel, sql, json, parquet,…). Importing data from each of these
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
data you need are available in Pandas.

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

Pandas provides plotting your data out of the box, using the power of Matplotlib. You can pick the plot type (scatter, bar, boxplot,...)
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

Pandas has great support for time series and has an extensive set of tools for working with dates, times, and time-indexed data.

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

Data sets do not only contain numerical data. Pandas provides a wide range of functions to cleaning textual data and extract useful information from it.

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

Currently working with other software for data manipulation in a tabular format? You're probably familiar to typical
data operations and know *what* to do with your tabular data, but lacking the syntax to execute these operations. Get to know
the pandas syntax by looking for equivalents from the software you already know:

.. container:: container

    .. container:: row

        .. container:: col-lg-6 col-md-6 col-sm-6 col-xs-12 d-flex

            .. container:: card text-center intro-card shadow

                .. image:: ../_static/logo_r.svg
                    :class: card-img-top
                    :height: 72px
                    :alt: R project logo

                .. container:: card-body flex-fill

                    The `R programming language <https://www.r-project.org/>`_ provides the ``data.frame`` data
                    structure and multiple packages, such as `tidyverse <https://www.tidyverse.org/>`_ use and
                    extend ``data.frame`` s for convenient data handling functionalities similar to pandas.

                    .. container:: custom-button

                        :ref:`Learn more <compare_with_r>`

        .. container:: col-lg-6 col-md-6 col-sm-6 col-xs-12 d-flex

            .. container:: card text-center intro-card shadow

                .. image:: ../_static/logo_sql.svg
                    :class: card-img-top
                    :height: 72px
                    :alt: SQL logo

                .. container:: card-body flex-fill

                    Already familiar to ``SELECT``, ``GROUP BY``, ``JOIN``,...?
                    Most of these SQL manipulations do have equivalents in pandas.

                    .. container:: custom-button

                        :ref:`Learn more <compare_with_sql>`

        .. container:: col-lg-6 col-md-6 col-sm-6 col-xs-12 d-flex

            .. container:: card text-center intro-card shadow

                .. image:: ../_static/logo_stata.svg
                    :class: card-img-top
                    :height: 52px
                    :alt: Stata logo

                .. container:: card-body flex-fill

                    The ``data set`` included in the `STATA <https://en.wikipedia.org/wiki/Stata>`_ statistical software
                    suite corresponds to the pandas ``data.frame``. Many of the operations known from STATA have an equivalent
                    in pandas.

                    .. container:: custom-button

                        :ref:`Learn more <compare_with_stata>`

        .. container:: col-lg-6 col-md-6 col-sm-6 col-xs-12 d-flex

            .. container:: card text-center intro-card shadow

                .. image:: ../_static/logo_sas.svg
                    :class: card-img-top
                    :height: 52px
                    :alt: SAS logo

                .. container:: card-body flex-fill

                    The  `SAS statistical software suite <https://en.wikipedia.org/wiki/SAS_(software)>`_
                    provides the ``data set`` corresponding to the pandas ``data.frame``. Also vectorized operations,
                    filtering, string processing operations,... from SAS have similar functions in pandas.

                    .. container:: custom-button

                        :ref:`Learn more <compare_with_sas>`

Community tutorials
-------------------

The community produces a wide variety of tutorials available online. Some of the
material is enlisted in the community contributed :ref:`tutorials`.


.. If you update this toctree, also update the manual toctree in the
   main index.rst.template

.. toctree::
    :maxdepth: 2
    :hidden:

    install
    overview
    10min
    intro_tutorials/index
    basics
    dsintro
    comparison/index
    tutorials
