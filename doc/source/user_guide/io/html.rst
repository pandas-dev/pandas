====
HTML
====

.. _io.read_html:

Reading HTML content
''''''''''''''''''''

.. warning::

   We **highly encourage** you to read the :ref:`HTML Table Parsing gotchas <io.html.gotchas>`
   below regarding the issues surrounding the BeautifulSoup4/html5lib/lxml parsers.

The top-level :func:`~pandas.io.html.read_html` function can accept an HTML
string/file/URL and will parse HTML tables into list of pandas ``DataFrames``.
Let's look at a few examples.

.. note::

   ``read_html`` returns a ``list`` of ``DataFrame`` objects, even if there is
   only a single table contained in the HTML content.

Read a URL with no options:

.. code-block:: ipython

   In [320]: url = "https://www.fdic.gov/resources/resolutions/bank-failures/failed-bank-list"

   In [321]: pd.read_html(url)
   Out[321]:
   [                         Bank NameBank           CityCity StateSt  ...              Acquiring InstitutionAI Closing DateClosing FundFund
    0                    Almena State Bank             Almena      KS  ...                          Equity Bank    October 23, 2020    10538
    1           First City Bank of Florida  Fort Walton Beach      FL  ...            United Fidelity Bank, fsb    October 16, 2020    10537
    2                 The First State Bank      Barboursville      WV  ...                       MVB Bank, Inc.       April 3, 2020    10536
    3                   Ericson State Bank            Ericson      NE  ...           Farmers and Merchants Bank   February 14, 2020    10535
    4     City National Bank of New Jersey             Newark      NJ  ...                      Industrial Bank    November 1, 2019    10534
    ..                                 ...                ...     ...  ...                                  ...                 ...      ...
    558                 Superior Bank, FSB           Hinsdale      IL  ...                Superior Federal, FSB       July 27, 2001     6004
    559                Malta National Bank              Malta      OH  ...                    North Valley Bank         May 3, 2001     4648
    560    First Alliance Bank & Trust Co.         Manchester      NH  ...  Southern New Hampshire Bank & Trust    February 2, 2001     4647
    561  National State Bank of Metropolis         Metropolis      IL  ...              Banterra Bank of Marion   December 14, 2000     4646
    562                   Bank of Honolulu           Honolulu      HI  ...                   Bank of the Orient    October 13, 2000     4645

    [563 rows x 7 columns]]

.. note::

   The data from the above URL changes every Monday so the resulting data above may be slightly different.

Read a URL while passing headers alongside the HTTP request:

.. code-block:: ipython

   In [322]: url = 'https://www.sump.org/notes/request/' # HTTP request reflector

   In [323]: pd.read_html(url)
   Out[323]:
   [                   0                    1
    0     Remote Socket:  51.15.105.256:51760
    1  Protocol Version:             HTTP/1.1
    2    Request Method:                  GET
    3       Request URI:      /notes/request/
    4     Request Query:                  NaN,
    0   Accept-Encoding:             identity
    1              Host:         www.sump.org
    2        User-Agent:    Python-urllib/3.8
    3        Connection:                close]

   In [324]: headers = {
      .....:    'User-Agent':'Mozilla Firefox v14.0',
      .....:    'Accept':'application/json',
      .....:    'Connection':'keep-alive',
      .....:    'Auth':'Bearer 2*/f3+fe68df*4'
      .....: }

   In [325]: pd.read_html(url, storage_options=headers)
   Out[325]:
   [                   0                    1
    0     Remote Socket:  51.15.105.256:51760
    1  Protocol Version:             HTTP/1.1
    2    Request Method:                  GET
    3       Request URI:      /notes/request/
    4     Request Query:                  NaN,
    0        User-Agent: Mozilla Firefox v14.0
    1    AcceptEncoding:   gzip,  deflate,  br
    2            Accept:      application/json
    3        Connection:             keep-alive
    4              Auth:  Bearer 2*/f3+fe68df*4]

.. note::

   We see above that the headers we passed are reflected in the HTTP request.

Read in the content of the file from the above URL and pass it to ``read_html``
as a string:

.. ipython:: python

   html_str = """
            <table>
                <tr>
                    <th>A</th>
                    <th colspan="1">B</th>
                    <th rowspan="1">C</th>
                </tr>
                <tr>
                    <td>a</td>
                    <td>b</td>
                    <td>c</td>
                </tr>
            </table>
        """

   with open("tmp.html", "w") as f:
       f.write(html_str)
   df = pd.read_html("tmp.html")
   df[0]

.. ipython:: python
   :suppress:

   import os
   os.remove("tmp.html")

You can even pass in an instance of ``StringIO`` if you so desire:

.. ipython:: python

   from io import StringIO

   dfs = pd.read_html(StringIO(html_str))
   dfs[0]

.. note::

   The following examples are not run by the IPython evaluator due to the fact
   that having so many network-accessing functions slows down the documentation
   build. If you spot an error or an example that doesn't run, please do not
   hesitate to report it over on `pandas GitHub issues page
   <https://github.com/pandas-dev/pandas/issues>`__.


Read a URL and match a table that contains specific text:

.. code-block:: python

   match = "Metcalf Bank"
   df_list = pd.read_html(url, match=match)

Specify a header row (by default ``<th>`` or ``<td>`` elements located within a
``<thead>`` are used to form the column index, if multiple rows are contained within
``<thead>`` then a MultiIndex is created); if specified, the header row is taken
from the data minus the parsed header elements (``<th>`` elements).

.. code-block:: python

   dfs = pd.read_html(url, header=0)

Specify an index column:

.. code-block:: python

   dfs = pd.read_html(url, index_col=0)

Specify a number of rows to skip:

.. code-block:: python

   dfs = pd.read_html(url, skiprows=0)

Specify a number of rows to skip using a list (``range`` works
as well):

.. code-block:: python

   dfs = pd.read_html(url, skiprows=range(2))

Specify an HTML attribute:

.. code-block:: python

   dfs1 = pd.read_html(url, attrs={"id": "table"})
   dfs2 = pd.read_html(url, attrs={"class": "sortable"})
   print(np.array_equal(dfs1[0], dfs2[0]))  # Should be True

Specify values that should be converted to NaN:

.. code-block:: python

   dfs = pd.read_html(url, na_values=["No Acquirer"])

Specify whether to keep the default set of NaN values:

.. code-block:: python

   dfs = pd.read_html(url, keep_default_na=False)

Specify converters for columns. This is useful for numerical text data that has
leading zeros.  By default columns that are numerical are cast to numeric
types and the leading zeros are lost. To avoid this, we can convert these
columns to strings.

.. code-block:: python

   url_mcc = "https://en.wikipedia.org/wiki/Mobile_country_code?oldid=899173761"
   dfs = pd.read_html(
       url_mcc,
       match="Telekom Albania",
       header=0,
       converters={"MNC": str},
   )

Use some combination of the above:

.. code-block:: python

   dfs = pd.read_html(url, match="Metcalf Bank", index_col=0)

Read in pandas ``to_html`` output (with some loss of floating point precision):

.. code-block:: python

   df = pd.DataFrame(np.random.randn(2, 2))
   s = df.to_html(float_format="{0:.40g}".format)
   dfin = pd.read_html(s, index_col=0)

The ``lxml`` backend will raise an error on a failed parse if that is the only
parser you provide. If you only have a single parser you can provide just a
string, but it is considered good practice to pass a list with one string if,
for example, the function expects a sequence of strings. You may use:

.. code-block:: python

   dfs = pd.read_html(url, "Metcalf Bank", index_col=0, flavor=["lxml"])

Or you could pass ``flavor='lxml'`` without a list:

.. code-block:: python

   dfs = pd.read_html(url, "Metcalf Bank", index_col=0, flavor="lxml")

However, if you have bs4 and html5lib installed and pass ``None`` or ``['lxml',
'bs4']`` then the parse will most likely succeed. Note that *as soon as a parse
succeeds, the function will return*.

.. code-block:: python

   dfs = pd.read_html(url, "Metcalf Bank", index_col=0, flavor=["lxml", "bs4"])

Links can be extracted from cells along with the text using ``extract_links="all"``.

.. ipython:: python

    from io import StringIO

    html_table = """
    <table>
      <tr>
        <th>GitHub</th>
      </tr>
      <tr>
        <td><a href="https://github.com/pandas-dev/pandas">pandas</a></td>
      </tr>
    </table>
    """

    df = pd.read_html(
        StringIO(html_table),
        extract_links="all"
    )[0]
    df
    df[("GitHub", None)]
    df[("GitHub", None)].str[1]

.. versionadded:: 1.5.0

.. _io.html:

Writing to HTML files
'''''''''''''''''''''

``DataFrame`` objects have an instance method ``to_html`` which renders the
contents of the ``DataFrame`` as an HTML table. The function arguments are as
in the method ``to_string`` described above.

.. note::

   Not all of the possible options for ``DataFrame.to_html`` are shown here for
   brevity's sake. See :func:`.DataFrame.to_html` for the
   full set of options.

.. note::

   In an HTML-rendering supported environment like a Jupyter Notebook, ``display(HTML(...))```
   will render the raw HTML into the environment.

.. ipython:: python

   from IPython.display import display, HTML

   df = pd.DataFrame(np.random.randn(2, 2))
   df
   html = df.to_html()
   print(html)  # raw html
   display(HTML(html))

The ``columns`` argument will limit the columns shown:

.. ipython:: python

   html = df.to_html(columns=[0])
   print(html)
   display(HTML(html))

``float_format`` takes a Python callable to control the precision of floating
point values:

.. ipython:: python

   html = df.to_html(float_format="{0:.10f}".format)
   print(html)
   display(HTML(html))


``bold_rows`` will make the row labels bold by default, but you can turn that
off:

.. ipython:: python

   html = df.to_html(bold_rows=False)
   print(html)
   display(HTML(html))


The ``classes`` argument provides the ability to give the resulting HTML
table CSS classes. Note that these classes are *appended* to the existing
``'dataframe'`` class.

.. ipython:: python

   print(df.to_html(classes=["awesome_table_class", "even_more_awesome_class"]))

The ``render_links`` argument provides the ability to add hyperlinks to cells
that contain URLs.

.. ipython:: python

   url_df = pd.DataFrame(
       {
           "name": ["Python", "pandas"],
           "url": ["https://www.python.org/", "https://pandas.pydata.org"],
       }
   )
   html = url_df.to_html(render_links=True)
   print(html)
   display(HTML(html))

Finally, the ``escape`` argument allows you to control whether the
"<", ">" and "&" characters escaped in the resulting HTML (by default it is
``True``). So to get the HTML without escaped characters pass ``escape=False``

.. ipython:: python

   df = pd.DataFrame({"a": list("&<>"), "b": np.random.randn(3)})

Escaped:

.. ipython:: python

   html = df.to_html()
   print(html)
   display(HTML(html))

Not escaped:

.. ipython:: python

   html = df.to_html(escape=False)
   print(html)
   display(HTML(html))

.. note::

   Some browsers may not show a difference in the rendering of the previous two
   HTML tables.


.. _io.html.gotchas:

HTML Table Parsing Gotchas
''''''''''''''''''''''''''

There are some versioning issues surrounding the libraries that are used to
parse HTML tables in the top-level pandas io function ``read_html``.

**Issues with** |lxml|_

* Benefits

    - |lxml|_ is very fast.

    - |lxml|_ requires Cython to install correctly.

* Drawbacks

    - |lxml|_ does *not* make any guarantees about the results of its parse
      *unless* it is given |svm|_.

    - In light of the above, we have chosen to allow you, the user, to use the
      |lxml|_ backend, but **this backend will use** |html5lib|_ if |lxml|_
      fails to parse

    - It is therefore *highly recommended* that you install both
      |BeautifulSoup4|_ and |html5lib|_, so that you will still get a valid
      result (provided everything else is valid) even if |lxml|_ fails.

**Issues with** |BeautifulSoup4|_ **using** |lxml|_ **as a backend**

* The above issues hold here as well since |BeautifulSoup4|_ is essentially
  just a wrapper around a parser backend.

**Issues with** |BeautifulSoup4|_ **using** |html5lib|_ **as a backend**

* Benefits

    - |html5lib|_ is far more lenient than |lxml|_ and consequently deals
      with *real-life markup* in a much saner way rather than just, e.g.,
      dropping an element without notifying you.

    - |html5lib|_ *generates valid HTML5 markup from invalid markup
      automatically*. This is extremely important for parsing HTML tables,
      since it guarantees a valid document. However, that does NOT mean that
      it is "correct", since the process of fixing markup does not have a
      single definition.

    - |html5lib|_ is pure Python and requires no additional build steps beyond
      its own installation.

* Drawbacks

    - The biggest drawback to using |html5lib|_ is that it is slow as
      molasses.  However consider the fact that many tables on the web are not
      big enough for the parsing algorithm runtime to matter. It is more
      likely that the bottleneck will be in the process of reading the raw
      text from the URL over the web, i.e., IO (input-output). For very large
      tables, this might not be true.


.. |svm| replace:: **strictly valid markup**
.. _svm: https://validator.w3.org/docs/help.html#validation_basics

.. |html5lib| replace:: **html5lib**
.. _html5lib: https://github.com/html5lib/html5lib-python

.. |BeautifulSoup4| replace:: **BeautifulSoup4**
.. _BeautifulSoup4: https://www.crummy.com/software/BeautifulSoup

.. |lxml| replace:: **lxml**
.. _lxml: https://lxml.de
