.. _io.other:

================================
Community-supported file formats
================================

pandas itself only supports IO with a limited set of file formats that map
cleanly to its tabular data model. For reading and writing other file formats
into and from pandas, we recommend these packages from the broader community.

.. _io.bigquery:

Google BigQuery
'''''''''''''''

The pandas-gbq_ package provides functionality to read/write from Google BigQuery.

.. _pandas-gbq: https://pandas-gbq.readthedocs.io/en/latest/

netCDF
''''''

xarray_ provides data structures inspired by the pandas ``DataFrame`` for working
with multi-dimensional datasets, with a focus on the netCDF file format and
easy conversion to and from pandas.

.. _xarray: https://xarray.pydata.org/en/stable/
