.. currentmodule:: pandas
.. _gotchas:

.. ipython:: python
   :suppress:

   import os
   import numpy as np
   from pandas import *
   randn = np.random.randn
   np.set_printoptions(precision=4, suppress=True)

*******************
Caveats and Gotchas
*******************

``NaN``, Integer ``NA`` values and ``NA`` type promotions
---------------------------------------------------------

Choice of ``NA`` representation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

For lack of ``NA`` (missing) support from the ground up in NumPy and Python in
general, we were given the difficult choice between either

- A *masked array* solution: an array of data and an array of boolean values
  indicating whether a value
- Using a special sentinel value, bit pattern, or set of sentinel values to
  denote ``NA`` across the dtypes

For many reasons we chose the latter. After years of production use it has
proven, at least in my opinion, to be the best decision given the state of
affairs in NumPy and Python in general. The special value ``NaN``
(Not-A-Number) is used everywhere as the ``NA`` value, and there are API
functions ``isnull`` and ``notnull`` which can be used across the dtypes to
detect NA values.

However, it comes with it a couple of trade-offs which I most certainly have
not ignored.

.. _gotchas.intna:

Support for integer ``NA``
~~~~~~~~~~~~~~~~~~~~~~~~~~

In the absence of high performance ``NA`` support being built into NumPy from
the ground up, the primary casualty is the ability to represent NAs in integer
arrays. For example:

.. ipython:: python

   s = Series([1, 2, 3, 4, 5], index=list('abcde'))
   s
   s.dtype

   s2 = s.reindex(['a', 'b', 'c', 'f', 'u'])
   s2
   s2.dtype

This trade-off is made largely for memory and performance reasons, and also so
that the resulting Series continues to be "numeric". One possibility is to use
``dtype=object`` arrays instead.

``NA`` type promotions
~~~~~~~~~~~~~~~~~~~~~~

When introducing NAs into an existing Series or DataFrame via ``reindex`` or
some other means, boolean and integer types will be promoted to a different
dtype in order to store the NAs. These are summarized by this table:

.. csv-table::
   :header: "Typeclass","Promotion dtype for storing NAs"
   :widths: 40,60

   ``floating``, no change
   ``object``, no change
   ``integer``, cast to ``float64``
   ``boolean``, cast to ``object``

While this may seem like a heavy trade-off, in practice I have found very few
cases where this is an issue in practice. Some explanation for the motivation
here in the next section.

Why not make NumPy like R?
~~~~~~~~~~~~~~~~~~~~~~~~~~

Many people have suggested that NumPy should simply emulate the ``NA`` support
present in the more domain-specific statistical programming langauge `R
<http://r-project.org>`__. Part of the reason is the NumPy type hierarchy:

.. csv-table::
   :header: "Typeclass","Dtypes"
   :widths: 30,70
   :delim: |

   ``numpy.floating`` | ``float16, float32, float64, float128``
   ``numpy.integer`` | ``int8, int16, int32, int64``
   ``numpy.unsignedinteger`` | ``uint8, uint16, uint32, uint64``
   ``numpy.object_`` | ``object_``
   ``numpy.bool_`` | ``bool_``
   ``numpy.character`` | ``string_, unicode_``

The R language, by contrast, only has a handful of built-in data types:
``integer``, ``numeric`` (floating-point), ``character``, and
``boolean``. ``NA`` types are implemented by reserving special bit patterns for
each type to be used as the missing value. While doing this with the full NumPy
type hierarchy would be possible, it would be a more substantial trade-off
(especially for the 8- and 16-bit data types) and implementation undertaking.

An alternate approach is that of using masked arrays. A masked array is an
array of data with an associated boolean *mask* denoting whether each value
should be considered ``NA`` or not. I am personally not in love with this
approach as I feel that overall it places a fairly heavy burden on the user and
the library implementer. Additionally, it exacts a fairly high performance cost
when working with numerical data compared with the simple approach of using
``NaN``. Thus, I have chosen the Pythonic "practicality beats purity" approach
and traded integer ``NA`` capability for a much simpler approach of using a
special value in float and object arrays to denote ``NA``, and promoting
integer arrays to floating when NAs must be introduced.

Integer indexing
----------------

Label-based indexing with integer axis labels is a thorny topic. It has been
discussed heavily on mailing lists and among various members of the scientific
Python community. In pandas, our general viewpoint is that labels matter more
than integer locations. Therefore, with an integer axis index *only*
label-based indexing is possible with the standard tools like ``.ix``. The
following code will generate exceptions:

.. code-block:: python

   s = Series(range(5))
   s[-1]
   df = DataFrame(np.random.randn(5, 4))
   df
   df.ix[-2:]

This deliberate decision was made to prevent ambiguities and subtle bugs (many
users reported finding bugs when the API change was made to stop "falling back"
on position-based indexing).

Label-based slicing conventions
-------------------------------

Non-monotonic indexes require exact matches
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Endpoints are inclusive
~~~~~~~~~~~~~~~~~~~~~~~

Compared with standard Python sequence slicing in which the slice endpoint is
not inclusive, label-based slicing in pandas **is inclusive**. The primary
reason for this is that it is often not possible to easily determine the
"successor" or next element after a particular label in an index. For example,
consider the following Series:

.. ipython:: python

   s = Series(randn(6), index=list('abcdef'))
   s

Suppose we wished to slice from ``c`` to ``e``, using integers this would be

.. ipython:: python

   s[2:5]

However, if you only had ``c`` and ``e``, determining the next element in the
index can be somewhat complicated. For example, the following does not work:

::

    s.ix['c':'e'+1]

A very common use case is to limit a time series to start and end at two
specific dates. To enable this, we made the design design to make label-based
slicing include both endpoints:

.. ipython:: python

    s.ix['c':'e']

This is most definitely a "practicality beats purity" sort of thing, but it is
something to watch out for if you expect label-based slicing to behave exactly
in the way that standard Python integer slicing works.

Miscellaneous indexing gotchas
------------------------------

Reindex versus ix gotchas
~~~~~~~~~~~~~~~~~~~~~~~~~

Many users will find themselves using the ``ix`` indexing capabilities as a
concise means of selecting data from a pandas object:

.. ipython:: python

   df = DataFrame(randn(6, 4), columns=['one', 'two', 'three', 'four'],
                  index=list('abcdef'))
   df
   df.ix[['b', 'c', 'e']]

This is, of course, completely equivalent *in this case* to using th
``reindex`` method:

.. ipython:: python

   df.reindex(['b', 'c', 'e'])

Some might conclude that ``ix`` and ``reindex`` are 100% equivalent based on
this. This is indeed true **except in the case of integer indexing**. For
example, the above operation could alternately have been expressed as:

.. ipython:: python

   df.ix[[1, 2, 4]]

If you pass ``[1, 2, 4]`` to ``reindex`` you will get another thing entirely:

.. ipython:: python

   df.reindex([1, 2, 4])

So it's important to remember that ``reindex`` is **strict label indexing
only**. This can lead to some potentially surprising results in pathological
cases where an index contains, say, both integers and strings:

.. ipython:: python

   s = Series([1, 2, 3], index=['a', 0, 1])
   s
   s.ix[[0, 1]]
   s.reindex([0, 1])

Because the index in this case does not contain solely integers, ``ix`` falls
back on integer indexing. By contrast, ``reindex`` only looks for the values
passed in the index, thus finding the integers ``0`` and ``1``. While it would
be possible to insert some logic to check whether a passed sequence is all
contained in the index, that logic would exact a very high cost in large data
sets.

Reindex potentially changes underlying Series dtype
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The use of ``reindex_like`` can potentially change the dtype of a ``Series``.

.. code-block:: python

   series = pandas.Series([1, 2, 3])
   x = pandas.Series([True])
   x.dtype
   x = pandas.Series([True]).reindex_like(series)
   x.dtype

This is because ``reindex_like`` silently inserts ``NaNs`` and the ``dtype``
changes accordingly.  This can cause some issues when using ``numpy`` ``ufuncs``
such as ``numpy.logical_and``.

See the `this old issue <https://github.com/pydata/pandas/issues/2388>`__ for a more
detailed discussion.

Timestamp limitations
---------------------

Minimum and maximum timestamps
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Since pandas represents timestamps in nanosecond resolution, the timespan that
can be represented using a 64-bit integer is limited to approximately 584 years:

.. ipython:: python

   begin = Timestamp.min
   begin

   end = Timestamp.max
   end

If you need to represent time series data outside the nanosecond timespan, use
PeriodIndex:

.. ipython:: python

   span = period_range('1215-01-01', '1381-01-01', freq='D')
   span


Parsing Dates from Text Files
-----------------------------

When parsing multiple text file columns into a single date column, the new date
column is prepended to the data and then `index_col` specification is indexed off
of the new set of columns rather than the original ones:

.. ipython:: python
   :suppress:

   data =  ("KORD,19990127, 19:00:00, 18:56:00, 0.8100\n"
            "KORD,19990127, 20:00:00, 19:56:00, 0.0100\n"
            "KORD,19990127, 21:00:00, 20:56:00, -0.5900\n"
            "KORD,19990127, 21:00:00, 21:18:00, -0.9900\n"
            "KORD,19990127, 22:00:00, 21:56:00, -0.5900\n"
            "KORD,19990127, 23:00:00, 22:56:00, -0.5900")

   with open('tmp.csv', 'w') as fh:
       fh.write(data)

.. ipython:: python

   print open('tmp.csv').read()

   date_spec = {'nominal': [1, 2], 'actual': [1, 3]}
   df = read_csv('tmp.csv', header=None,
                 parse_dates=date_spec,
                 keep_date_col=True,
                 index_col=0)

   # index_col=0 refers to the combined column "nominal" and not the original
   # first column of 'KORD' strings

   df

.. ipython:: python
   :suppress:

   os.remove('tmp.csv')


Differences with NumPy
----------------------
For Series and DataFrame objects, ``var`` normalizes by ``N-1`` to produce
unbiased estimates of the sample variance, while NumPy's ``var`` normalizes
by N, which measures the variance of the sample. Note that ``cov``
normalizes by ``N-1`` in both pandas and NumPy.


Thread-safety
-------------

As of pandas 0.11, pandas is not 100% thread safe. The known issues relate to
the ``DataFrame.copy`` method. If you are doing a lot of copying of DataFrame
objects shared among threads, we recommend holding locks inside the threads
where the data copying occurs.

See `this link <http://stackoverflow.com/questions/13592618/python-pandas-dataframe-thread-safe>`__
for more information.

.. _html-gotchas:

HTML Table Parsing
------------------
There are some versioning issues surrounding the libraries that are used to
parse HTML tables in the top-level pandas io function ``read_html``.

**Issues with** |lxml|_

   * Benefits

     * |lxml|_ is very fast

     * |lxml|_ requires Cython to install correctly.

   * Drawbacks

     * |lxml|_ does *not* make any guarantees about the results of it's parse
       *unless* it is given |svm|_.

     * In light of the above, we have chosen to allow you, the user, to use the
       |lxml|_ backend, but **this backend will use** |html5lib|_ if |lxml|_
       fails to parse

     * It is therefore *highly recommended* that you install both
       |BeautifulSoup4|_ and |html5lib|_, so that you will still get a valid
       result (provided everything else is valid) even if |lxml|_ fails.

**Issues with** |BeautifulSoup4|_ **using** |lxml|_ **as a backend**

   * The above issues hold here as well since |BeautifulSoup4|_ is essentially
     just a wrapper around a parser backend.

**Issues with** |BeautifulSoup4|_ **using** |html5lib|_ **as a backend**

   * Benefits

     * |html5lib|_ is far more lenient than |lxml|_ and consequently deals
       with *real-life markup* in a much saner way rather than just, e.g.,
       dropping an element without notifying you.

     * |html5lib|_ *generates valid HTML5 markup from invalid markup
       automatically*. This is extremely important for parsing HTML tables,
       since it guarantees a valid document. However, that does NOT mean that
       it is "correct", since the process of fixing markup does not have a
       single definition.

     * |html5lib|_ is pure Python and requires no additional build steps beyond
       its own installation.

   * Drawbacks

     * The biggest drawback to using |html5lib|_ is that it is slow as
       molasses.  However consider the fact that many tables on the web are not
       big enough for the parsing algorithm runtime to matter. It is more
       likely that the bottleneck will be in the process of reading the raw
       text from the url over the web, i.e., IO (input-output). For very large
       tables, this might not be true.

**Issues with using** |Anaconda|_

   * `Anaconda`_ ships with `lxml`_ version 3.2.0; the following workaround for
     `Anaconda`_ was successfully used to deal with the versioning issues
     surrounding `lxml`_ and `BeautifulSoup4`_.

   .. note::

      Unless you have *both*:

         * A strong restriction on the upper bound of the runtime of some code
           that incorporates :func:`~pandas.io.html.read_html`
         * Complete knowledge that the HTML you will be parsing will be 100%
           valid at all times

      then you should install `html5lib`_ and things will work swimmingly
      without you having to muck around with `conda`. If you want the best of
      both worlds then install both `html5lib`_ and `lxml`_. If you do install
      `lxml`_ then you need to perform the following commands to ensure that
      lxml will work correctly:

      .. code-block:: sh
         
         # remove the included version
         conda remove lxml

         # install the latest version of lxml
         pip install 'git+git://github.com/lxml/lxml.git'

         # install the latest version of beautifulsoup4
         pip install 'bzr+lp:beautifulsoup'

      Note that you need `bzr <http://bazaar.canonical.com/en>`_ and `git
      <http://git-scm.com>`_ installed to perform the last two operations.

.. |svm| replace:: **strictly valid markup**
.. _svm: http://validator.w3.org/docs/help.html#validation_basics

.. |html5lib| replace:: **html5lib**
.. _html5lib: https://github.com/html5lib/html5lib-python

.. |BeautifulSoup4| replace:: **BeautifulSoup4**
.. _BeautifulSoup4: http://www.crummy.com/software/BeautifulSoup

.. |lxml| replace:: **lxml**
.. _lxml: http://lxml.de

.. |Anaconda| replace:: **Anaconda**
.. _Anaconda: https://store.continuum.io/cshop/anaconda


Byte-Ordering Issues
--------------------
Occasionally you may have to deal with data that were created on a machine with
a different byte order than the one on which you are running Python. To deal
with this issue you should convert the underlying NumPy array to the native
system byte order *before* passing it to Series/DataFrame/Panel constructors
using something similar to the following:

.. ipython:: python

   x = np.array(range(10), '>i4') # big endian
   newx = x.byteswap().newbyteorder() # force native byteorder
   s = Series(newx)

See `the NumPy documentation on byte order
<http://docs.scipy.org/doc/numpy/user/basics.byteswapping.html>`__ for more
details.
