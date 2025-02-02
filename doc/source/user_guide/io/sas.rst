.. _io.sas:

.. _io.sas_reader:

===========
SAS formats
===========

The top-level function :func:`read_sas` can read (but not write) SAS
XPORT (.xpt) and SAS7BDAT (.sas7bdat) format files.

SAS files only contain two value types: ASCII text and floating point
values (usually 8 bytes but sometimes truncated).  For xport files,
there is no automatic type conversion to integers, dates, or
categoricals.  For SAS7BDAT files, the format codes may allow date
variables to be automatically converted to dates.  By default the
whole file is read and returned as a ``DataFrame``.

Specify a ``chunksize`` or use ``iterator=True`` to obtain reader
objects (``XportReader`` or ``SAS7BDATReader``) for incrementally
reading the file.  The reader objects also have attributes that
contain additional information about the file and its variables.

Read a SAS7BDAT file:

.. code-block:: python

    df = pd.read_sas("sas_data.sas7bdat")

Obtain an iterator and read an XPORT file 100,000 lines at a time:

.. code-block:: python

    def do_something(chunk):
        pass


    with pd.read_sas("sas_xport.xpt", chunk=100000) as rdr:
        for chunk in rdr:
            do_something(chunk)

The specification_ for the xport file format is available from the SAS
web site.

.. _specification: https://support.sas.com/content/dam/SAS/support/en/technical-papers/record-layout-of-a-sas-version-5-or-6-data-set-in-sas-transport-xport-format.pdf

No official documentation is available for the SAS7BDAT format.
