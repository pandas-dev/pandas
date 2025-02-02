.. _io:

===============================
IO tools (text, CSV, HDF5, ...)
===============================

.. toctree::
    :maxdepth: 1
    :hidden:

    csv
    json
    html
    latex
    xml
    clipboard
    excel
    hdf5
    feather
    parquet
    orc
    stata
    sas
    spss
    pickling
    sql
    community_packages

The pandas I/O API is a set of top level ``reader`` functions accessed like
:func:`pandas.read_csv` that generally return a pandas object. The corresponding
``writer`` functions are object methods that are accessed like
:meth:`DataFrame.to_csv`. Below is a table containing available ``readers`` and
``writers``.

.. csv-table::
    :header: "Format Type", "Data Description", "Reader", "Writer"
    :widths: 30, 100, 60, 60

    text,`CSV <https://en.wikipedia.org/wiki/Comma-separated_values>`__, :ref:`read_csv<io.read_csv_table>`, :ref:`to_csv<io.store_in_csv>`
    text,Fixed-Width Text File, :ref:`read_fwf<io.fwf_reader>` , NA
    text,`JSON <https://www.json.org/>`__, :ref:`read_json<io.json_reader>`, :ref:`to_json<io.json_writer>`
    text,`HTML <https://en.wikipedia.org/wiki/HTML>`__, :ref:`read_html<io.read_html>`, :ref:`to_html<io.html>`
    text,`LaTeX <https://en.wikipedia.org/wiki/LaTeX>`__, :ref:`Styler.to_latex<io.latex>` , NA
    text,`XML <https://www.w3.org/standards/xml/core>`__, :ref:`read_xml<io.read_xml>`, :ref:`to_xml<io.xml>`
    text, Local clipboard, :ref:`read_clipboard<io.clipboard>`, :ref:`to_clipboard<io.clipboard>`
    binary,`MS Excel <https://en.wikipedia.org/wiki/Microsoft_Excel>`__ , :ref:`read_excel<io.excel_reader>`, :ref:`to_excel<io.excel_writer>`
    binary,`OpenDocument <http://opendocumentformat.org>`__, :ref:`read_excel<io.ods>`, NA
    binary,`HDF5 Format <https://support.hdfgroup.org/HDF5/whatishdf5.html>`__, :ref:`read_hdf<io.hdf5>`, :ref:`to_hdf<io.hdf5>`
    binary,`Feather Format <https://github.com/wesm/feather>`__, :ref:`read_feather<io.feather>`, :ref:`to_feather<io.feather>`
    binary,`Parquet Format <https://parquet.apache.org/>`__, :ref:`read_parquet<io.parquet>`, :ref:`to_parquet<io.parquet>`
    binary,`ORC Format <https://orc.apache.org/>`__, :ref:`read_orc<io.orc>`, :ref:`to_orc<io.orc>`
    binary,`Stata <https://en.wikipedia.org/wiki/Stata>`__, :ref:`read_stata<io.stata_reader>`, :ref:`to_stata<io.stata_writer>`
    binary,`SAS <https://en.wikipedia.org/wiki/SAS_(software)>`__, :ref:`read_sas<io.sas_reader>` , NA
    binary,`SPSS <https://en.wikipedia.org/wiki/SPSS>`__, :ref:`read_spss<io.spss_reader>` , NA
    binary,`Python Pickle Format <https://docs.python.org/3/library/pickle.html>`__, :ref:`read_pickle<io.pickle>`, :ref:`to_pickle<io.pickle>`
    SQL,`SQL <https://en.wikipedia.org/wiki/SQL>`__, :ref:`read_sql<io.sql>`,:ref:`to_sql<io.sql>`

See also this list of :ref:`community-supported packages<io.other>` offering support for other file formats.


.. _io.perf:

Performance considerations
--------------------------

This is an informal comparison of various IO methods, using pandas
0.24.2. Timings are machine dependent and small differences should be
ignored.

.. code-block:: ipython

   In [1]: sz = 1000000
   In [2]: df = pd.DataFrame({'A': np.random.randn(sz), 'B': [1] * sz})

   In [3]: df.info()
   <class 'pandas.DataFrame'>
   RangeIndex: 1000000 entries, 0 to 999999
   Data columns (total 2 columns):
   A    1000000 non-null float64
   B    1000000 non-null int64
   dtypes: float64(1), int64(1)
   memory usage: 15.3 MB

The following test functions will be used below to compare the performance of several IO methods:

.. code-block:: python



   import numpy as np

   import os

   sz = 1000000
   df = pd.DataFrame({"A": np.random.randn(sz), "B": [1] * sz})

   sz = 1000000
   np.random.seed(42)
   df = pd.DataFrame({"A": np.random.randn(sz), "B": [1] * sz})


   def test_sql_write(df):
       if os.path.exists("test.sql"):
           os.remove("test.sql")
       sql_db = sqlite3.connect("test.sql")
       df.to_sql(name="test_table", con=sql_db)
       sql_db.close()


   def test_sql_read():
       sql_db = sqlite3.connect("test.sql")
       pd.read_sql_query("select * from test_table", sql_db)
       sql_db.close()


   def test_hdf_fixed_write(df):
       df.to_hdf("test_fixed.hdf", key="test", mode="w")


   def test_hdf_fixed_read():
       pd.read_hdf("test_fixed.hdf", "test")


   def test_hdf_fixed_write_compress(df):
       df.to_hdf("test_fixed_compress.hdf", key="test", mode="w", complib="blosc")


   def test_hdf_fixed_read_compress():
       pd.read_hdf("test_fixed_compress.hdf", "test")


   def test_hdf_table_write(df):
       df.to_hdf("test_table.hdf", key="test", mode="w", format="table")


   def test_hdf_table_read():
       pd.read_hdf("test_table.hdf", "test")


   def test_hdf_table_write_compress(df):
       df.to_hdf(
           "test_table_compress.hdf", key="test", mode="w", complib="blosc", format="table"
       )


   def test_hdf_table_read_compress():
       pd.read_hdf("test_table_compress.hdf", "test")


   def test_csv_write(df):
       df.to_csv("test.csv", mode="w")


   def test_csv_read():
       pd.read_csv("test.csv", index_col=0)


   def test_feather_write(df):
       df.to_feather("test.feather")


   def test_feather_read():
       pd.read_feather("test.feather")


   def test_pickle_write(df):
       df.to_pickle("test.pkl")


   def test_pickle_read():
       pd.read_pickle("test.pkl")


   def test_pickle_write_compress(df):
       df.to_pickle("test.pkl.compress", compression="xz")


   def test_pickle_read_compress():
       pd.read_pickle("test.pkl.compress", compression="xz")


   def test_parquet_write(df):
       df.to_parquet("test.parquet")


   def test_parquet_read():
       pd.read_parquet("test.parquet")

When writing, the top three functions in terms of speed are ``test_feather_write``, ``test_hdf_fixed_write`` and ``test_hdf_fixed_write_compress``.

.. code-block:: ipython

   In [4]: %timeit test_sql_write(df)
   3.29 s ± 43.2 ms per loop (mean ± std. dev. of 7 runs, 1 loop each)

   In [5]: %timeit test_hdf_fixed_write(df)
   19.4 ms ± 560 µs per loop (mean ± std. dev. of 7 runs, 1 loop each)

   In [6]: %timeit test_hdf_fixed_write_compress(df)
   19.6 ms ± 308 µs per loop (mean ± std. dev. of 7 runs, 10 loops each)

   In [7]: %timeit test_hdf_table_write(df)
   449 ms ± 5.61 ms per loop (mean ± std. dev. of 7 runs, 1 loop each)

   In [8]: %timeit test_hdf_table_write_compress(df)
   448 ms ± 11.9 ms per loop (mean ± std. dev. of 7 runs, 1 loop each)

   In [9]: %timeit test_csv_write(df)
   3.66 s ± 26.2 ms per loop (mean ± std. dev. of 7 runs, 1 loop each)

   In [10]: %timeit test_feather_write(df)
   9.75 ms ± 117 µs per loop (mean ± std. dev. of 7 runs, 100 loops each)

   In [11]: %timeit test_pickle_write(df)
   30.1 ms ± 229 µs per loop (mean ± std. dev. of 7 runs, 10 loops each)

   In [12]: %timeit test_pickle_write_compress(df)
   4.29 s ± 15.9 ms per loop (mean ± std. dev. of 7 runs, 1 loop each)

   In [13]: %timeit test_parquet_write(df)
   67.6 ms ± 706 µs per loop (mean ± std. dev. of 7 runs, 10 loops each)

When reading, the top three functions in terms of speed are ``test_feather_read``, ``test_pickle_read`` and
``test_hdf_fixed_read``.


.. code-block:: ipython

   In [14]: %timeit test_sql_read()
   1.77 s ± 17.7 ms per loop (mean ± std. dev. of 7 runs, 1 loop each)

   In [15]: %timeit test_hdf_fixed_read()
   19.4 ms ± 436 µs per loop (mean ± std. dev. of 7 runs, 10 loops each)

   In [16]: %timeit test_hdf_fixed_read_compress()
   19.5 ms ± 222 µs per loop (mean ± std. dev. of 7 runs, 10 loops each)

   In [17]: %timeit test_hdf_table_read()
   38.6 ms ± 857 µs per loop (mean ± std. dev. of 7 runs, 10 loops each)

   In [18]: %timeit test_hdf_table_read_compress()
   38.8 ms ± 1.49 ms per loop (mean ± std. dev. of 7 runs, 10 loops each)

   In [19]: %timeit test_csv_read()
   452 ms ± 9.04 ms per loop (mean ± std. dev. of 7 runs, 1 loop each)

   In [20]: %timeit test_feather_read()
   12.4 ms ± 99.7 µs per loop (mean ± std. dev. of 7 runs, 100 loops each)

   In [21]: %timeit test_pickle_read()
   18.4 ms ± 191 µs per loop (mean ± std. dev. of 7 runs, 100 loops each)

   In [22]: %timeit test_pickle_read_compress()
   915 ms ± 7.48 ms per loop (mean ± std. dev. of 7 runs, 1 loop each)

   In [23]: %timeit test_parquet_read()
   24.4 ms ± 146 µs per loop (mean ± std. dev. of 7 runs, 10 loops each)


The files ``test.pkl.compress``, ``test.parquet`` and ``test.feather`` took the least space on disk (in bytes).

.. code-block:: none

    29519500 Oct 10 06:45 test.csv
    16000248 Oct 10 06:45 test.feather
    8281983  Oct 10 06:49 test.parquet
    16000857 Oct 10 06:47 test.pkl
    7552144  Oct 10 06:48 test.pkl.compress
    34816000 Oct 10 06:42 test.sql
    24009288 Oct 10 06:43 test_fixed.hdf
    24009288 Oct 10 06:43 test_fixed_compress.hdf
    24458940 Oct 10 06:44 test_table.hdf
    24458940 Oct 10 06:44 test_table_compress.hdf
