{{ header }}

.. _getting_started:

===============
Getting started
===============

Installation
------------

Before you can use pandas, you’ll need to get it installed.

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

Pandas is part of the `Anaconda <http://docs.continuum.io/anaconda/>`__ distribution and can be 
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

Pandas can be installed via pip from `PyPI <https://pypi.org/project/pandas>`__.

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

TODO


.. _comingfrom:

Coming from...
--------------

Currently working with other software for data manipulation in a tabular format? You're probably familiar to typical
data operations and know *what* to do with your tabular data, but lacking the syntax to execute these operations. Get to know
the pandas syntax by looking for equivalents from the software you already know:

.. raw:: html

    <div class="container">
        <div class="row">
            <div class="col-lg-6 col-md-6 col-sm-6 col-xs-12 d-flex">
                <div class="card text-center intro-card shadow">
                <img src="_static/logo_r.svg" class="card-img-top" alt="R project logo" height="72"> 
                <div class="card-body flex-fill">
                    <p class="card-text">The <a href="https://www.r-project.org/">R programming language</a> provides the <code>data.frame</code> data structure and multiple packages, 
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
                <img src="_static/logo_sql.svg" class="card-img-top" alt="SQL logo" height="72"> 
                <div class="card-body flex-fill">
                    <p class="card-text">Already familiar to <code>SELECT</code>, <code>GROUP BY</code>, <code>JOIN</code>,...? 
                    Most of these SQL manipulations do have equivalents in pandas.</p>

.. container:: custom-button

    :ref:`Learn more <compare_with_sql>`

.. raw:: html

                    </div>
                    </div>
                </div>
                <div class="col-lg-6 col-md-6 col-sm-6 col-xs-12 d-flex">
                    <div class="card text-center intro-card shadow">
                    <img src="_static/logo_stata.svg" class="card-img-top" alt="STATA logo" height="52"> 
                    <div class="card-body flex-fill">
                        <p class="card-text">The <code>data set</code> included in the 
                            <a href="https://en.wikipedia.org/wiki/Stata">STATA</a> statistical software suite corresponds 
                            to the pandas <code>data.frame</code>. Many of the operations known from STATA have an equivalent
                            in pandas.</p>

.. container:: custom-button
    
    :ref:`Learn more <compare_with_stata>`

.. raw:: html

                    </div>
                    </div>
                </div>
                <div class="col-lg-6 col-md-6 col-sm-6 col-xs-12 d-flex">
                    <div class="card text-center intro-card shadow">
                    <img src="_static/logo_sas.svg" class="card-img-top" alt="SAS logo" height="52"> 
                    <div class="card-body flex-fill">
                        <p class="card-text">The  <a href="https://en.wikipedia.org/wiki/SAS_(software)">SAS</a> statistical software suite 
                            also provides the <code>data set</code> corresponding to the pandas <code>data.frame</code>. 
                            Also vectorized operations, filtering, string processing operations,... from SAS have similar 
                            functions in pandas.</p>

.. container:: custom-button

    :ref:`Learn more <compare_with_sas>`

.. raw:: html

                    </div>
                    </div>
                </div>           
        </div>		
    </div>

Community tutorials
-------------------

The community produces a wide variety of tutorials available online. Some of the 
material is enlisted in the community contributed :ref:`tutorials`. 


.. If you update this toctree, also update the manual toctree in the
   main index.rst.template

.. toctree::
    :maxdepth: 2

    install
    overview
    10min
    intro_tutorials/index
    basics
    dsintro
    comparison/index
    tutorials
