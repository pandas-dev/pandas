.. _contributing.docs:

Contributing to the documentation
=================================

If you're not the developer type, contributing to the documentation is still
of huge value. You don't even have to be an expert on
*pandas* to do so! Something as simple as rewriting small passages for clarity
as you reference the docs is a simple but effective way to contribute. The
next person to read that passage will be in your debt!

Actually, there are sections of the docs that are worse off by being written
by experts. If something in the docs doesn't make sense to you, updating the
relevant section after you figure it out is a simple way to ensure it will
help the next person.

.. contents:: Table of contents:
   :local:


About the pandas documentation
------------------------------

The documentation is written in **reStructuredText**, which is almost like writing
in plain English, and built using `Sphinx <http://sphinx.pocoo.org/>`__. The
Sphinx Documentation has an excellent `introduction to reST
<http://sphinx.pocoo.org/rest.html>`__. Review the Sphinx docs to perform more
complex changes to the documentation as well.

Some other important things to know about the docs:

- The pandas documentation consists of two parts: the docstrings in the code
  itself and the docs in this folder ``pandas/doc/``.

  The docstrings provide a clear explanation of the usage of the individual
  functions, while the documentation in this folder consists of tutorial-like
  overviews per topic together with some other information (what's new,
  installation, etc).

- The docstrings follow the **Numpy Docstring Standard** which is used widely
  in the Scientific Python community. This standard specifies the format of
  the different sections of the docstring. See `this document
  <https://github.com/numpy/numpy/blob/master/doc/HOWTO_DOCUMENT.rst.txt>`_
  for a detailed explanation, or look at some of the existing functions to
  extend it in a similar manner.

- The tutorials make heavy use of the `ipython directive
  <http://matplotlib.org/sampledoc/ipython_directive.html>`_ sphinx extension.
  This directive lets you put code in the documentation which will be run
  during the doc build. For example:

  ::

      .. ipython:: python

          x = 2
          x**3

  will be rendered as

  ::

      In [1]: x = 2

      In [2]: x**3
      Out[2]: 8

  This means that almost all code examples in the docs are always run (and the
  output saved) during the doc build. This way, they will always be up to date,
  but it makes the doc building a bit more complex.


How to build the pandas documentation
-------------------------------------

Requirements
^^^^^^^^^^^^

To build the pandas docs there are some extra requirements: you will need to
have ``sphinx`` and ``ipython`` installed. `numpydoc
<https://github.com/numpy/numpydoc>`_ is used to parse the docstrings that
follow the Numpy Docstring Standard (see above), but you don't need to install
this because a local copy of ``numpydoc`` is included in the pandas source
code.

Furthermore, it is recommended to have all `optional dependencies
<http://pandas.pydata.org/pandas-docs/dev/install.html#optional-dependencies>`_
installed. This is not needed, but be aware that you will see some error
messages. Because all the code in the documentation is executed during the doc
build, the examples using this optional dependencies will generate errors.
Run ``pd.show_versions()`` to get an overview of the installed version of all
dependencies.

.. warning::

   Sphinx version >= 1.2.2 or the older 1.1.3 is required.

Building pandas
^^^^^^^^^^^^^^^

For a step-by-step overview on how to set up your environment, to work with
the pandas code and git, see `the developer pages
<http://pandas.pydata.org/developers.html#working-with-the-code>`_.
When you start to work on some docs, be sure to update your code to the latest
development version ('master')::

    git fetch upstream
    git rebase upstream/master

Often it will be necessary to rebuild the C extension after updating::

    python setup.py build_ext --inplace

Building the documentation
^^^^^^^^^^^^^^^^^^^^^^^^^^

So how do you build the docs? Navigate to your local folder
``pandas/doc/`` directory in the console and run::

    python make.py html

And then you can find the html output in the folder ``pandas/doc/build/html/``.

The first time it will take quite a while, because it has to run all the code
examples in the documentation and build all generated docstring pages.
In subsequent evocations, sphinx will try to only build the pages that have
been modified.

If you want to do a full clean build, do::

    python make.py clean
    python make.py build


Starting with 0.13.1 you can tell ``make.py`` to compile only a single section
of the docs, greatly reducing the turn-around time for checking your changes.
You will be prompted to delete `.rst` files that aren't required, since the
last committed version can always be restored from git.

::

    #omit autosummary and API section
    python make.py clean
    python make.py --no-api

    # compile the docs with only a single
    # section, that which is in indexing.rst
    python make.py clean
    python make.py --single indexing

For comparison, a full doc build may take 10 minutes. a ``-no-api`` build
may take 3 minutes and a single section may take 15 seconds.

Where to start?
---------------

There are a number of issues listed under `Docs
<https://github.com/pandas-dev/pandas/issues?labels=Docs&sort=updated&state=open>`_
and `Good as first PR
<https://github.com/pandas-dev/pandas/issues?labels=Good+as+first+PR&sort=updated&state=open>`_
where you could start out.

Or maybe you have an idea of your own, by using pandas, looking for something
in the documentation and thinking 'this can be improved', let's do something
about that!

Feel free to ask questions on `mailing list
<https://groups.google.com/forum/?fromgroups#!forum/pydata>`_ or submit an
issue on Github.
