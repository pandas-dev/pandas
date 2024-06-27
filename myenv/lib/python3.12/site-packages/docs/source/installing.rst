Installing airspeed velocity
============================

**airspeed velocity** is known to work on Linux, Mac OS-X, and Windows.
It is known to work with Python 3.7 and higher.
It works also with PyPy.

**airspeed velocity** is a standard Python package, and the latest
released version may be `installed in the standard
way from PyPI <https://packaging.python.org/tutorials/installing-packages/>`__::

    pip install asv

The development version can be installed by cloning the source
repository and running ``pip install .`` inside it, or by ``pip
install git+https://github.com/airspeed-velocity/asv``.

The requirements should be automatically installed.  If they aren't
installed automatically, for example due to networking restrictions,
the ``python`` requirements are as noted in the ``pyproject.toml``.

For managing the environments, one of the following packages is required:

- `libmambapy <https://mamba.readthedocs.io/en/latest/python_api.html>`__,
  which is typically part of ``mamba``. In this case ``conda`` must be present too.

- `virtualenv <http://virtualenv.org/>`__, which is required since
  venv is not compatible with other versions of Python.

- An `anaconda <https://store.continuum.io/cshop/anaconda/>`__ or
  `miniconda <http://conda.pydata.org/miniconda.html>`__
  installation, with the ``conda`` command available on your path.

.. note::

   ``libmambapy`` is the fastest for situations where non-pythonic
   dependencies are required. Anaconda or miniconda is slower but
   still preferred if the project involves a lot of compiled C/C++
   extensions and are available in the ``conda`` repository, since
   ``conda`` will be able to fetch precompiled binaries for these
   dependencies in many cases. Using ``virtualenv``, dependencies
   without precompiled wheels usually have to be compiled every
   time the environments are set up.

Optional optimization
---------------------

If your project being benchmarked contains C, C++, Objective-C or
Cython, consider installing ``ccache``.  `ccache
<https://ccache.samba.org/>`__ is a compiler cache that speeds up
compilation time when the same objects are repeatedly compiled.  In
**airspeed velocity**, the project being benchmarked is recompiled at
many different points in its history, often with only minor changes to
the source code, so ``ccache`` can help speed up the total benchmarking
time considerably.

Running the self-tests
----------------------

The self tests are based on `pytest <https://docs.pytest.org/>`__.  If you
don't have it installed, and you have a connection to the Internet, it
will be installed automatically.

To run **airspeed velocity**'s self tests::

    pytest
