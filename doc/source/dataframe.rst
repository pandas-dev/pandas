.. _dataframe:

.. currentmodule:: pandas

************************
DataFrame and DataMatrix
************************

**DataFrame** is a both 2D-matrix- and dict-like object implementing a
named collection of identically-indexed Series objects. There is
another important class called **DataMatrix** which has almost
identical behavior to DataFrame, differing only in the internal
implementation. Rather than create confusion, we will focus on
DataFrame for this guide and explain later why to use one over the
other depending on your application.

For many users, DataFrame will be the most commonly used object in
pandas, as it can serve as the primary container for a full data set
of interest. As we will see below, it enables operations involving
multiple time series or cross sections with ease.

.. class:: DataFrame

   :Parameters:
       **data** : dict or 2-D ndarray
           * If the dict contains Series, an index need not be specified.
             The resulting frame will be the union of all the contained
             Series indices.

       **index** : {array_like}
           Explicit index to conform to, required for ndarray data argument

       **columns** : {array_like}
           Explicit set of columns to include, required for ndarray data argument

       **dtype** : Python type alias or :class:`~numpy.dtype`
           Type to attempt to cast data to

:class:`~pandas.DataMatrix` has a similar constructor and can be used interchangeably.


Basics
------

.. note::

    Unlike Series, neither DataFrame nor DataMatrix is a subclass of
    numpy.ndarray.

The canonical DataFrame containing time series data takes this form,
which will be used for many examples to follow:

::

    >>> from pandas import *; from numpy.random import randn

    >>> index = DateRange('1/1/2009', '12/1/2009', timeRule='EOM')
    >>> N = len(index)
    >>> data = {
	'A' : randn(N),
	'B' : randn(N),
	'C' : randn(N),
    }

    >>> df = DataFrame(data, index=index)

    >>> print df

			   A              B              C
    2009-01-30 00:00:00    -0.367014      1.29942        1.40773
    2009-02-27 00:00:00    0.347326       0.651661       0.143376
    2009-03-31 00:00:00    -0.677813      0.40488        -0.463113
    2009-04-30 00:00:00    0.125062       1.09505        -1.03278
    2009-05-29 00:00:00    0.979307       0.149356       0.708128
    2009-06-30 00:00:00    -1.24432       0.420788       -1.01815
    2009-07-31 00:00:00    0.536261       -0.276357      -0.227469
    2009-08-31 00:00:00    0.0603968      -1.42112       -0.271767
    2009-09-30 00:00:00    0.537841       -0.361833      -0.0488729
    2009-10-30 00:00:00    -0.335999      2.54742        -0.878263
    2009-11-30 00:00:00    -0.568216      -0.557347      -1.58623

The **info** method provides a summary of a DataFrame object and will
be printed by default when the frame is very large:

::

    >>> df.info()
    Index: 11 entries, 2009-01-30 00:00:00 to 2009-11-30 00:00:00
    Columns:
    A    11  non-null values
    B    11  non-null values
    C    11  non-null values

The DataFrame's index and columns can be accessed by the **index**
attribute and **cols** method, respectively:

::

    >>> df.index
    Index([2009-01-30 00:00:00, 2009-02-27 00:00:00, 2009-03-31 00:00:00,
	   2009-04-30 00:00:00, 2009-05-29 00:00:00, 2009-06-30 00:00:00,
	   2009-07-31 00:00:00, 2009-08-31 00:00:00, 2009-09-30 00:00:00,
	   2009-10-30 00:00:00, 2009-11-30 00:00:00], dtype=object)

    >>> df.cols()
    ['A', 'B', 'C']


.. autosummary::
   :toctree: generated/

   DataFrame.cols
   DataFrame.info

Construction
------------

There are many ways to create a DataFrame:

   * From a dict of ndarrays or Series
   * From a 2D ndarray plus corresponding row and column labels
   * From a NumPy structured (record) array
   * From a nested dictionary

.. autosummary::
   :toctree: generated/

   DataFrame.__init__
   DataFrame.fromRecords

Indexing
--------

.. note::

    I have eschewed "cute" indexing schemes in the interest of being
    explicit and maintaining DataFrame's status as a "dict of
    Series". However, it is always desirable to keep the interface
    simple and intuitive.


DataFrame's basic indexing accesses the **columns** by name, producing
Series:

::

    >>> df['A']
    2009-01-30 00:00:00    -0.367013536107
    2009-02-27 00:00:00    0.347325830717
    2009-03-31 00:00:00    -0.677812757268
    2009-04-30 00:00:00    0.125061634713
    2009-05-29 00:00:00    0.979307492892
    2009-06-30 00:00:00    -1.2443243316
    2009-07-31 00:00:00    0.536260924391
    2009-08-31 00:00:00    0.060396849998
    2009-09-30 00:00:00    0.53784064627
    2009-10-30 00:00:00    -0.335999254912
    2009-11-30 00:00:00    -0.568216482894


If you add a Series to the frame, it will be automatically conformed
to the frame's index:

::

    >>> df['D'] = df['A'][:5]
    >>> df
    			   A              B              C              D
    2009-01-30 00:00:00    -0.367014      1.29942        1.40773        -0.367014
    2009-02-27 00:00:00    0.347326       0.651661       0.143376       0.347326
    2009-03-31 00:00:00    -0.677813      0.40488        -0.463113      -0.677813
    2009-04-30 00:00:00    0.125062       1.09505        -1.03278       0.125062
    2009-05-29 00:00:00    0.979307       0.149356       0.708128       0.979307
    2009-06-30 00:00:00    -1.24432       0.420788       -1.01815       nan
    2009-07-31 00:00:00    0.536261       -0.276357      -0.227469      nan
    2009-08-31 00:00:00    0.0603968      -1.42112       -0.271767      nan
    2009-09-30 00:00:00    0.537841       -0.361833      -0.0488729     nan
    2009-10-30 00:00:00    -0.335999      2.54742        -0.878263      nan
    2009-11-30 00:00:00    -0.568216      -0.557347      -1.58623       nan

Columns can be deleted or popped as with a dict:

::

    >>> del df['C']
    >>> B = df.pop('B')
    >>> df.info()
    Index: 11 entries, 2009-01-30 00:00:00 to 2009-11-30 00:00:00
    Columns:
    A    11  non-null values
    D    5  non-null values


New items in the DataFrame do not need to already be Series. They can
also be an ndarray of the right length or a scalar value:

::

    >>> df['N'] = np.arange(len(df))
    >>> df['S'] = 5
			   A              D              N              S
    2009-01-30 00:00:00    -0.367014      -0.367014      0              5
    2009-02-27 00:00:00    0.347326       0.347326       1              5
    2009-03-31 00:00:00    -0.677813      -0.677813      2              5
    2009-04-30 00:00:00    0.125062       0.125062       3              5
    2009-05-29 00:00:00    0.979307       0.979307       4              5
    2009-06-30 00:00:00    -1.24432       nan            5              5
    2009-07-31 00:00:00    0.536261       nan            6              5
    2009-08-31 00:00:00    0.0603968      nan            7              5
    2009-09-30 00:00:00    0.537841       nan            8              5
    2009-10-30 00:00:00    -0.335999      nan            9              5
    2009-11-30 00:00:00    -0.568216      nan            10             5

To be consistent with this dict-like interface, the *__contains__*
method considers the columns:

::

    >>> 'A' in df
    True

.. autosummary::
   :toctree: generated/

   DataFrame.__contains__
   DataFrame.__getitem__
   DataFrame.__delitem__
   DataFrame.pop

Retrieving cross sections, transposing
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

It is often desirable to retrieve all data associated with a
particular index values (we have been calling this a *cross
section*). Rather than use *__getitem__* and extra notation to do
this, DataFrame has the **xs** method:

::

    >>> df.xs(datetime(2009, 8, 31))
    A    0.060396849998
    D    nan
    N    7
    S    5

If the cross sections in a data set are of the most interest, it is
also possible to transpose a DataFrame using the same syntax as an
ndarray:

::

    >>> dftrans = df.T
    >>> dftrans
    <class 'pandas.core.frame.DataFrame'>
    Index: 4 entries, A to S
    Columns:
    2009-01-30 00:00:00    4  non-null values
    2009-02-27 00:00:00    4  non-null values
    2009-03-31 00:00:00    4  non-null values
    2009-04-30 00:00:00    4  non-null values
    2009-05-29 00:00:00    4  non-null values
    2009-06-30 00:00:00    3  non-null values
    2009-07-31 00:00:00    3  non-null values
    2009-08-31 00:00:00    3  non-null values
    2009-09-30 00:00:00    3  non-null values
    2009-10-30 00:00:00    3  non-null values
    2009-11-30 00:00:00    3  non-null values

    >>> dftrans[datetime(2009, 9, 30)]
    A    0.53784064627
    D    nan
    N    8
    S    5

Slicing ranges
~~~~~~~~~~~~~~

Similar to Python lists and ndarrays, for convenience DataFrame
supports slicing:

::

    >>> df[:2]
			   A              D              N              S
    2009-01-30 00:00:00    -0.367014      -0.367014      0              5
    2009-02-27 00:00:00    0.347326       0.347326       1              5

    >>> df[::-1]
			   A              D              N              S
    2009-11-30 00:00:00    -0.568216      nan            10             5
    2009-10-30 00:00:00    -0.335999      nan            9              5
    2009-09-30 00:00:00    0.537841       nan            8              5
    2009-08-31 00:00:00    0.0603968      nan            7              5
    2009-07-31 00:00:00    0.536261       nan            6              5
    2009-06-30 00:00:00    -1.24432       nan            5              5
    2009-05-29 00:00:00    0.979307       0.979307       4              5
    2009-04-30 00:00:00    0.125062       0.125062       3              5
    2009-03-31 00:00:00    -0.677813      -0.677813      2              5
    2009-02-27 00:00:00    0.347326       0.347326       1              5
    2009-01-30 00:00:00    -0.367014      -0.367014      0              5

    >>> df[-3:].T
	 2009-09-30     2009-10-30     2009-11-30
    A    0.537841       -0.335999      -0.568216
    D    nan            nan            nan
    N    8              9              10
    S    5              5              5

I do not recommend making heavy use of this functionality but rather
using it as a convenience for interactive programming (useful for
seeing the "head" or "tail" of a large DataFrame as in the last
example).

Boolean
~~~~~~~

As another indexing convenience, it is possible to use boolean
indexing to select rows of a DataFrame:

::

    >>> df[df['A'] > 0.5]
			   A              D              N              S
    2009-05-29 00:00:00    0.979307       0.979307       4              5
    2009-07-31 00:00:00    0.536261       nan            6              5
    2009-09-30 00:00:00    0.537841       nan            8              5

As we will see later on, the same operation could be accomplished by
reindexing. However, the syntax would be more verbose; hence, the
inclusion of this indexing method.

Arithmetic
----------

.. seealso:: :ref:`Series arithmetic <series.arithmetic>`

Binary operations with DataFrame have similar index-matching behavior
as with Series objects. The addition is the matching of column names
between DataFrame objects or Series. We will detail how the
interactions work in each case.

DataFrame and DataFrame
~~~~~~~~~~~~~~~~~~~~~~~

When combining two DataFrames, both index and column values must match
for two values to be combined. If there is no match for a particular
(index, column) pair, the result for that location will be NaN. To
illustrate, let's return to a similar example from the beginning of
the tutorial:

::

    >>> df
			   A              B              C
    2009-01-30 00:00:00    -0.173487      -0.330054      2.45767
    2009-02-27 00:00:00    -1.70517       -1.34422       0.45781
    2009-03-31 00:00:00    0.517951       0.437294       0.625021
    2009-04-30 00:00:00    1.13914        0.976763       0.871074
    2009-05-29 00:00:00    -0.263249      -1.55445       0.386744
    2009-06-30 00:00:00    0.994217       -0.15012       -0.444482
    2009-07-31 00:00:00    1.51264        -1.13902       0.846015
    2009-08-31 00:00:00    0.323804       -0.793455      -1.97154
    2009-09-30 00:00:00    -0.0450052     0.404083       0.588554
    2009-10-30 00:00:00    0.268981       -0.20756       -0.328061
    2009-11-30 00:00:00    0.471714       -0.0450022     -0.280202

    >>> df + df[:7]
			   A              B              C
    2009-01-30 00:00:00    -0.346974      -0.660108      4.91534
    2009-02-27 00:00:00    -3.41035       -2.68845       0.915621
    2009-03-31 00:00:00    1.0359         0.874588       1.25004
    2009-04-30 00:00:00    2.27828        1.95353        1.74215
    2009-05-29 00:00:00    -0.526498      -3.10889       0.773488
    2009-06-30 00:00:00    1.98843        -0.30024       -0.888965
    2009-07-31 00:00:00    3.02528        -2.27803       1.69203
    2009-08-31 00:00:00    nan            nan            nan
    2009-09-30 00:00:00    nan            nan            nan
    2009-10-30 00:00:00    nan            nan            nan
    2009-11-30 00:00:00    nan            nan            nan

In this first example, we can see that the indices have been combined
together, and the portion where dates are missing in one of the frames
has resulted in all NaN values. The resulting columns will also be the
union of the frames' columns:

::

    >>> df2 = df.copy()
    >>> df2['D'] = 5
    >>> del df2['A']

    >>> df + df2[::2]
			   A              B              C              D
    2009-01-30 00:00:00    nan            -0.660108      4.91534        nan
    2009-02-27 00:00:00    nan            nan            nan            nan
    2009-03-31 00:00:00    nan            0.874588       1.25004        nan
    2009-04-30 00:00:00    nan            nan            nan            nan
    2009-05-29 00:00:00    nan            -3.10889       0.773488       nan
    2009-06-30 00:00:00    nan            nan            nan            nan
    2009-07-31 00:00:00    nan            -2.27803       1.69203        nan
    2009-08-31 00:00:00    nan            nan            nan            nan
    2009-09-30 00:00:00    nan            0.808167       1.17711        nan
    2009-10-30 00:00:00    nan            nan            nan            nan
    2009-11-30 00:00:00    nan            -0.0900043     -0.560404      nan

Here, neither **A** nor **D** was in both frames: they appear in the
result but are all NaN. An argument could be made to exclude these
columns, but it is very frequently meaningful to know that there was
no overlap in a particular portion of the data set.

DataFrame and Series
~~~~~~~~~~~~~~~~~~~~

The choice of behavior between DataFrame and Series was somewhat
arbitrary. Since the **columns** of the DataFrame are viewed as its
keys, and the **index** values of a Series are viewed as *its* keys,
the default behavior is to match the frame columns on the series
index.

::

    >>> df - df.xs(df.index[5])
			   A              B              C
    2009-01-30 00:00:00    -1.1677        -0.179934      2.90215
    2009-02-27 00:00:00    -2.69939       -1.1941        0.902293
    2009-03-31 00:00:00    -0.476266      0.587414       1.0695
    2009-04-30 00:00:00    0.144924       1.12688        1.31556
    2009-05-29 00:00:00    -1.25747       -1.40433       0.831226
    2009-06-30 00:00:00    0              0              0
    2009-07-31 00:00:00    0.518423       -0.988895      1.2905
    2009-08-31 00:00:00    -0.670413      -0.643335      -1.52705
    2009-09-30 00:00:00    -1.03922       0.554203       1.03304
    2009-10-30 00:00:00    -0.725236      -0.0574399     0.116421
    2009-11-30 00:00:00    -0.522503      0.105118       0.16428

However, the user very frequently will want to subtract (or add,
divide, multiply, ...) a TimeSeries from a DataFrame representing a
collection of TimeSeries. Since this is so common, the DataFrame will
inspect the input Series and its own index to see if this is what the
user intended.

::

    >>> df - df['A']
			   A              B              C
    2009-01-30 00:00:00    0              -0.156567      2.63116
    2009-02-27 00:00:00    0              0.360951       2.16298
    2009-03-31 00:00:00    0              -0.0806571     0.10707
    2009-04-30 00:00:00    0              -0.162378      -0.268067
    2009-05-29 00:00:00    0              -1.2912        0.649993
    2009-06-30 00:00:00    0              -1.14434       -1.4387
    2009-07-31 00:00:00    0              -2.65166       -0.666625
    2009-08-31 00:00:00    0              -1.11726       -2.29534
    2009-09-30 00:00:00    0              0.449089       0.633559
    2009-10-30 00:00:00    0              -0.476541      -0.597042
    2009-11-30 00:00:00    0              -0.516716      -0.751916

Note that the same result could have been obtained by writing:

::

    >>> (df.T - df['A']).T

but this is fairly awkward (and relatively high cost due to two
transpose operations).

DataFrame and scalar value
~~~~~~~~~~~~~~~~~~~~~~~~~~

Scalar operations work just as expected:

::

    >>> df * 5 + 2
			   A              B              C
    2009-01-30 00:00:00    1.13256        0.349729       14.2884
    2009-02-27 00:00:00    -6.52587       -4.72111       4.28905
    2009-03-31 00:00:00    4.58976        4.18647        5.12511
    2009-04-30 00:00:00    7.69571        6.88382        6.35537
    2009-05-29 00:00:00    0.683756       -5.77223       3.93372
    2009-06-30 00:00:00    6.97108        1.2494         -0.222412
    2009-07-31 00:00:00    9.5632         -3.69508       6.23008
    2009-08-31 00:00:00    3.61902        -1.96728       -7.85768
    2009-09-30 00:00:00    1.77497        4.02042        4.94277
    2009-10-30 00:00:00    3.34491        0.962201       0.359695
    2009-11-30 00:00:00    4.35857        1.77499        0.59899

    >>> 1 / df
    ...
    >>> df ** 4
    ...

Basic statistical functions
---------------------------

.. seealso:: :ref:`Series statistical methods <series.statistics>`

Being used to working with ndarrays, we would like to have the same
sorts of basic descriptive statistics implemented for
DataFrame. Additionally, and similar to NumPy MaskedArray objects,
these should be able to handle missing data.

An analogous set of methods to those in Series are provided to compute
common moments and other aggregate statistics. These come with the
addition that the aggregation can be over either axis of the
DataFrame. By default the statistic will be computed for each column
(axis 0):

::

    >>> df
			   A              B              C              D
    2009-01-30 00:00:00    -0.173487      -0.330054      2.45767        -0.173487
    2009-02-27 00:00:00    -1.70517       -1.34422       0.45781        -1.70517
    2009-03-31 00:00:00    0.517951       0.437294       0.625021       0.517951
    2009-04-30 00:00:00    1.13914        0.976763       0.871074       1.13914
    2009-05-29 00:00:00    -0.263249      -1.55445       0.386744       -0.263249
    2009-06-30 00:00:00    0.994217       -0.15012       -0.444482      nan
    2009-07-31 00:00:00    1.51264        -1.13902       0.846015       nan
    2009-08-31 00:00:00    0.323804       -0.793455      -1.97154       nan
    2009-09-30 00:00:00    -0.0450052     0.404083       0.588554       nan
    2009-10-30 00:00:00    0.268981       -0.20756       -0.328061      nan
    2009-11-30 00:00:00    0.471714       -0.0450022     -0.280202      nan

    >>> df.mean()
    A	0.27650301895
    B	-0.340521353823
    C	0.291691664865
    D	-0.0969634747505

    >>> df.mean(axis=1)
    2009-01-30 00:00:00	0.445160910434
    2009-02-27 00:00:00	-1.07419006997
    2009-03-31 00:00:00	0.524554413383
    2009-04-30 00:00:00	1.03152984415
    2009-05-29 00:00:00	-0.42354987732
    2009-06-30 00:00:00	0.133204894302
    2009-07-31 00:00:00	0.406546724596
    2009-08-31 00:00:00	-0.813729414124
    2009-09-30 00:00:00	0.315877288659
    2009-10-30 00:00:00	-0.088879865383
    2009-11-30 00:00:00	0.048836496438

The other methods listed function similarly. Combining these methods
with the arithmetic functionality, we can very easily do things like
computing the cross-sectional or time series z-score:

::

    >>> (df - df.mean(1)) / df.std(1)    # cross-sectional
			   A              B              C              D
    2009-01-30 00:00:00    -0.839662      0.539768       1.13956        -0.839662
    2009-02-27 00:00:00    -0.615664      1.47701        -0.245682      -0.615664
    ...

    >>> (df - df.mean()) / df.std()      # time series
			   A              B              C              D
    2009-01-30 00:00:00    -1.79314       -0.173077      0.156054       -0.96237
    2009-02-27 00:00:00    -1.16502       1.7612         -1.15894       -0.437405
    ...

.. autosummary::
   :toctree: generated/

   DataFrame.count
   DataFrame.sum
   DataFrame.cumsum
   DataFrame.product
   DataFrame.mean
   DataFrame.median
   DataFrame.min
   DataFrame.max
   DataFrame.mad
   DataFrame.var
   DataFrame.std
   DataFrame.skew

Correlation
~~~~~~~~~~~

One related method computes the pairwise correlation of the columns,
taking care to compute the variances over the intersection of data:

::

    >>> df.corr()
	 A              B              C              D
    A    1              0.423766       -0.0985818     1
    B    0.423766       1              0.134803       0.839771
    C    -0.0985818     0.134803       1              0.129427
    D    1              0.839771       0.129427       1

Obviously other estimators of pairwise relationships could be computed
in the same way.

.. autosummary::
   :toctree: generated/

   DataFrame.corr

Function application
--------------------

You will often want to perform some other computation with the
DataFrame other than the statistical operators above. Likewise, you
may want to transform the data in some way (like taking the square
root or natural logarithm). To solve these problems, you should use
the **apply** function. In short, **apply** will call a function that
you pass on each row or column of the DataFrame and, depending on the
return type of the function, return a Series or DataFrame.

::

    >>> df.apply(np.log)
			   A              B              C
    2009-01-30 00:00:00    nan            nan            -0.573389
    2009-02-27 00:00:00    nan            -0.0123289     nan
    2009-03-31 00:00:00    0.0324552      -0.982897      -1.85985
    2009-04-30 00:00:00    nan            nan            nan
    2009-05-29 00:00:00    nan            nan            -0.70761
    2009-06-30 00:00:00    -1.01969       -0.598632      1.01358
    2009-07-31 00:00:00    -2.73689       nan            -0.051959
    2009-08-31 00:00:00    -1.90946       nan            -1.59065
    2009-09-30 00:00:00    -0.193634      -0.87138       -0.103805
    2009-10-30 00:00:00    -1.50367       -2.91441       -0.139465
    2009-11-30 00:00:00    -2.08095       nan            nan

    >>> df.apply(lambda x: np.sort(x)[-5:].mean())
    A	0.517625676559
    B	0.47682898773
    C	1.2079285542

    >>> df.apply(np.sum, axis=1)
    2009-01-30 00:00:00    -1.35501374903
    2009-02-27 00:00:00    -1.05605622659
    2009-03-31 00:00:00    1.56290932477
    2009-04-30 00:00:00    -2.1465835771
    2009-05-29 00:00:00    -1.14076063095
    2009-06-30 00:00:00    3.66570960365
    2009-07-31 00:00:00    0.328066060924
    2009-08-31 00:00:00    -0.130326883911
    2009-09-30 00:00:00    2.14373419389
    2009-10-30 00:00:00    1.14637153907
    2009-11-30 00:00:00    -1.12947529696

**apply** combined with some cleverness can be used to answer many
questions about a data set. For example, suppose we wanted to extract
the date where the maximum value for each column occurred:

::

    >>> df.apply(lambda x: df.index[x.valid().argmax()])
    A    2009-03-31 00:00:00
    B    2009-02-27 00:00:00
    C    2009-06-30 00:00:00

Another useful feature is the ability to pass Series methods to carry
out some Series operation on each columns or row:

::

    >>> df.apply(Series.interpolate)
    ...

applymap
~~~~~~~~

Since we don't always have vectorized functions to **apply**,
DataFrame has the method **applymap** which takes an elementwise
function. Obviously this will be fairly slow but can be useful:

::

    >>> df.applymap(lambda x: x if x > 0 else 0)
			   A              B              C              D
    2009-01-30 00:00:00    0              0              0.563612       0
    2009-02-27 00:00:00    0              0.987747       0              0
    2009-03-31 00:00:00    1.03299        0.374225       0.155696       1.03299
    2009-04-30 00:00:00    0              0              0              0
    2009-05-29 00:00:00    0              0              0.492821       0
    2009-06-30 00:00:00    0.360708       0.549563       2.75544        0
    2009-07-31 00:00:00    0.0647715      0              0.949368       0
    2009-08-31 00:00:00    0.148161       0              0.203793       0
    2009-09-30 00:00:00    0.82396        0.418374       0.901401       0
    2009-10-30 00:00:00    0.222313       0.0542359      0.869823       0
    2009-11-30 00:00:00    0.124812       0              0              0

Of course you could have accomplished this using a vectorized NumPy
function:

::

    >>> df.apply(lambda x: np.where(x > 0, x, 0))
    ...

Iterating
---------

DataFrame's iterator functions are consistent with the dict interface:

::

    >>> for column in df:
    ....    print column
    A
    C
    B
    D

    >>> for column, series in df.iteritems():
    ...

.. autosummary::
   :toctree: generated/

   DataFrame.cols
   DataFrame.iteritems

Reindexing
----------

.. seealso:: :ref:`Series reindexing <series.reindexing>`

Similar to Series, the **reindex** method conforms a DataFrame to a
new index or list of columns.

.. autosummary::
   :toctree: generated/

   DataFrame.reindex
   DataFrame.dropEmptyRows
   DataFrame.dropIncompleteRows
   DataFrame.merge
   DataFrame.fill
   DataFrame.filter

Sorting
-------

TODO

.. autosummary::
   :toctree: generated/

   DataFrame.sort

Converting to ndarray
---------------------

TODO

Joining / merging DataFrames
----------------------------

TODO

.. autosummary::
   :toctree: generated/

   DataFrame.join

TimeSeries-oriented methods
---------------------------

TODO

.. autosummary::
   :toctree: generated/

   DataFrame.asfreq
   DataFrame.truncate
   DataFrame.diff
   DataFrame.shift

Sorting
-------

TODO

.. autosummary::
   :toctree: generated/

   Series.argsort
   Series.sort
   Series.order

GroupBy functionality
---------------------

.. seealso:: :ref:`Series GroupBy <series.groupby>`

TODO

.. autosummary::
   :toctree: generated/

   DataFrame.groupby
   DataFrame.tgroupby

IO
--

TODO

.. autosummary::
   :toctree: generated/

   DataFrame.toCSV
   DataFrame.toString
   DataFrame.info

Miscellaneous
-------------

TODO

.. autosummary::
   :toctree: generated/

   DataFrame.append
   DataFrame.asMatrix
   DataFrame.values
   DataFrame.copy
   DataFrame.pivot
   DataFrame.T
   DataFrame.apply
   DataFrame.tapply
   DataFrame.applymap
   DataFrame.sort
   DataFrame.combineFirst
   DataFrame.combineAdd
   DataFrame.combineMult
   DataFrame.plot

DataFrame vs. DataMatrix
------------------------

DataFrame and DataMatrix differ only in their internal
implementation. In a nutshell:

  * DataFrame: dict of Series
  * DataMatrix: 2D homogeneous ndarray of values

