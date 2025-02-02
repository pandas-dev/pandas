.. _io.spss:

.. _io.spss_reader:

============
SPSS formats
============

The top-level function :func:`read_spss` can read (but not write) SPSS
SAV (.sav) and  ZSAV (.zsav) format files.

SPSS files contain column names. By default the
whole file is read, categorical columns are converted into ``pd.Categorical``,
and a ``DataFrame`` with all columns is returned.

Specify the ``usecols`` parameter to obtain a subset of columns. Specify ``convert_categoricals=False``
to avoid converting categorical columns into ``pd.Categorical``.

Read an SPSS file:

.. code-block:: python

    df = pd.read_spss("spss_data.sav")

Extract a subset of columns contained in ``usecols`` from an SPSS file and
avoid converting categorical columns into ``pd.Categorical``:

.. code-block:: python

    df = pd.read_spss(
        "spss_data.sav",
        usecols=["foo", "bar"],
        convert_categoricals=False,
    )

More information about the SAV and ZSAV file formats is available here_.

.. _here: https://www.ibm.com/docs/en/spss-statistics/22.0.0
