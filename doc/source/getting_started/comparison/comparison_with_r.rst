.. _compare_with_r:

{{ header }}

Comparison with R / R libraries
*******************************

Since pandas aims to provide a lot of the data manipulation and analysis
functionality that people use `R <https://www.r-project.org/>`__ for, this page
was started to provide a more detailed look at the `R language
<https://en.wikipedia.org/wiki/R_(programming_language)>`__ and its many third
party libraries as they relate to pandas. In comparisons with R and CRAN
libraries, we care about the following things:

* **Functionality / flexibility**: what can/cannot be done with each tool
* **Performance**: how fast are operations. Hard numbers/benchmarks are
  preferable
* **Ease-of-use**: Is one tool easier/harder to use (you may have to be
  the judge of this, given side-by-side code comparisons)

This page is also here to offer a bit of a translation guide for users of these
R packages.

For transfer of ``DataFrame`` objects from pandas to R, one option is to
use HDF5 files, see :ref:`io.external_compatibility` for an
example.


Quick reference
---------------

We'll start off with a quick reference guide pairing some common R
operations using `dplyr
<https://cran.r-project.org/package=dplyr>`__ with
pandas equivalents.


Querying, filtering, sampling
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

===========================================  ===========================================
R                                            pandas
===========================================  ===========================================
``dim(df)``                                  ``df.shape``
``head(df)``                                 ``df.head()``
``slice(df, 1:10)``                          ``df.iloc[:9]``
``filter(df, col1 == 1, col2 == 1)``         ``df.query('col1 == 1 & col2 == 1')``
``df[df$col1 == 1 & df$col2 == 1,]``         ``df[(df.col1 == 1) & (df.col2 == 1)]``
``select(df, col1, col2)``                   ``df[['col1', 'col2']]``
``select(df, col1:col3)``                    ``df.loc[:, 'col1':'col3']``
``select(df, -(col1:col3))``                 ``df.drop(cols_to_drop, axis=1)`` but see [#select_range]_
``distinct(select(df, col1))``               ``df[['col1']].drop_duplicates()``
``distinct(select(df, col1, col2))``         ``df[['col1', 'col2']].drop_duplicates()``
``sample_n(df, 10)``                         ``df.sample(n=10)``
``sample_frac(df, 0.01)``                    ``df.sample(frac=0.01)``
===========================================  ===========================================

.. [#select_range] R's shorthand for a subrange of columns
                   (``select(df, col1:col3)``) can be approached
                   cleanly in pandas, if you have the list of columns,
                   for example ``df[cols[1:3]]`` or
                   ``df.drop(cols[1:3])``, but doing this by column
                   name is a bit messy.


Sorting
~~~~~~~

===========================================  ===========================================
R                                            pandas
===========================================  ===========================================
``arrange(df, col1, col2)``                  ``df.sort_values(['col1', 'col2'])``
``arrange(df, desc(col1))``                  ``df.sort_values('col1', ascending=False)``
===========================================  ===========================================

Transforming
~~~~~~~~~~~~

===========================================  ===========================================
R                                            pandas
===========================================  ===========================================
``select(df, col_one = col1)``               ``df.rename(columns={'col1': 'col_one'})['col_one']``
``rename(df, col_one = col1)``               ``df.rename(columns={'col1': 'col_one'})``
``mutate(df, c=a-b)``                        ``df.assign(c=df['a']-df['b'])``
===========================================  ===========================================


Grouping and summarizing
~~~~~~~~~~~~~~~~~~~~~~~~

==============================================  ===========================================
R                                               pandas
==============================================  ===========================================
``summary(df)``                                 ``df.describe()``
``gdf <- group_by(df, col1)``                   ``gdf = df.groupby('col1')``
``summarise(gdf, avg=mean(col1, na.rm=TRUE))``  ``df.groupby('col1').agg({'col1': 'mean'})``
``summarise(gdf, total=sum(col1))``             ``df.groupby('col1').sum()``
==============================================  ===========================================


Base R
------

Slicing with R's |c|_
~~~~~~~~~~~~~~~~~~~~~

R makes it easy to access ``data.frame`` columns by name

.. code-block:: r

   df <- data.frame(a=rnorm(5), b=rnorm(5), c=rnorm(5), d=rnorm(5), e=rnorm(5))
   df[, c("a", "c", "e")]

or by integer location

.. code-block:: r

   df <- data.frame(matrix(rnorm(1000), ncol=100))
   df[, c(1:10, 25:30, 40, 50:100)]

Selecting multiple columns by name in pandas is straightforward

.. ipython:: python

   df = pd.DataFrame(np.random.randn(10, 3), columns=list("abc"))
   df[["a", "c"]]
   df.loc[:, ["a", "c"]]

Selecting multiple noncontiguous columns by integer location can be achieved
with a combination of the ``iloc`` indexer attribute and ``numpy.r_``.

.. ipython:: python

   named = list("abcdefg")
   n = 30
   columns = named + np.arange(len(named), n).tolist()
   df = pd.DataFrame(np.random.randn(n, n), columns=columns)

   df.iloc[:, np.r_[:10, 24:30]]

|aggregate|_
~~~~~~~~~~~~

In R you may want to split data into subsets and compute the mean for each.
Using a data.frame called ``df`` and splitting it into groups ``by1`` and
``by2``:

.. code-block:: r

   df <- data.frame(
     v1 = c(1,3,5,7,8,3,5,NA,4,5,7,9),
     v2 = c(11,33,55,77,88,33,55,NA,44,55,77,99),
     by1 = c("red", "blue", 1, 2, NA, "big", 1, 2, "red", 1, NA, 12),
     by2 = c("wet", "dry", 99, 95, NA, "damp", 95, 99, "red", 99, NA, NA))
   aggregate(x=df[, c("v1", "v2")], by=list(mydf2$by1, mydf2$by2), FUN = mean)

The :meth:`~pandas.DataFrame.groupby` method is similar to base R ``aggregate``
function.

.. ipython:: python

   df = pd.DataFrame(
       {
           "v1": [1, 3, 5, 7, 8, 3, 5, np.nan, 4, 5, 7, 9],
           "v2": [11, 33, 55, 77, 88, 33, 55, np.nan, 44, 55, 77, 99],
           "by1": ["red", "blue", 1, 2, np.nan, "big", 1, 2, "red", 1, np.nan, 12],
           "by2": [
               "wet",
               "dry",
               99,
               95,
               np.nan,
               "damp",
               95,
               99,
               "red",
               99,
               np.nan,
               np.nan,
           ],
       }
   )

   g = df.groupby(["by1", "by2"])
   g[["v1", "v2"]].mean()

For more details and examples see :ref:`the groupby documentation
<groupby.split>`.

|match|_
~~~~~~~~~~~~

A common way to select data in R is using ``%in%`` which is defined using the
function ``match``. The operator ``%in%`` is used to return a logical vector
indicating if there is a match or not:

.. code-block:: r

   s <- 0:4
   s %in% c(2,4)

The :meth:`~pandas.DataFrame.isin` method is similar to R ``%in%`` operator:

.. ipython:: python

   s = pd.Series(np.arange(5), dtype=np.float32)
   s.isin([2, 4])

The ``match`` function returns a vector of the positions of matches
of its first argument in its second:

.. code-block:: r

   s <- 0:4
   match(s, c(2,4))

For more details and examples see :ref:`the reshaping documentation
<indexing.basics.indexing_isin>`.

|tapply|_
~~~~~~~~~

``tapply`` is similar to ``aggregate``, but data can be in a ragged array,
since the subclass sizes are possibly irregular. Using a data.frame called
``baseball``, and retrieving information based on the array ``team``:

.. code-block:: r

   baseball <-
     data.frame(team = gl(5, 5,
                labels = paste("Team", LETTERS[1:5])),
                player = sample(letters, 25),
                batting.average = runif(25, .200, .400))

   tapply(baseball$batting.average, baseball.example$team,
          max)

In pandas we may use :meth:`~pandas.pivot_table` method to handle this:

.. ipython:: python

   import random
   import string

   baseball = pd.DataFrame(
       {
           "team": ["team %d" % (x + 1) for x in range(5)] * 5,
           "player": random.sample(list(string.ascii_lowercase), 25),
           "batting avg": np.random.uniform(0.200, 0.400, 25),
       }
   )

   baseball.pivot_table(values="batting avg", columns="team", aggfunc=np.max)

For more details and examples see :ref:`the reshaping documentation
<reshaping.pivot>`.

|subset|_
~~~~~~~~~~

The :meth:`~pandas.DataFrame.query` method is similar to the base R ``subset``
function. In R you might want to get the rows of a ``data.frame`` where one
column's values are less than another column's values:

.. code-block:: r

   df <- data.frame(a=rnorm(10), b=rnorm(10))
   subset(df, a <= b)
   df[df$a <= df$b,]  # note the comma

In pandas, there are a few ways to perform subsetting. You can use
:meth:`~pandas.DataFrame.query` or pass an expression as if it were an
index/slice as well as standard boolean indexing:

.. ipython:: python

   df = pd.DataFrame({"a": np.random.randn(10), "b": np.random.randn(10)})
   df.query("a <= b")
   df[df["a"] <= df["b"]]
   df.loc[df["a"] <= df["b"]]

For more details and examples see :ref:`the query documentation
<indexing.query>`.


|with|_
~~~~~~~~

An expression using a data.frame called ``df`` in R with the columns ``a`` and
``b`` would be evaluated using ``with`` like so:

.. code-block:: r

   df <- data.frame(a=rnorm(10), b=rnorm(10))
   with(df, a + b)
   df$a + df$b  # same as the previous expression

In pandas the equivalent expression, using the
:meth:`~pandas.DataFrame.eval` method, would be:

.. ipython:: python

   df = pd.DataFrame({"a": np.random.randn(10), "b": np.random.randn(10)})
   df.eval("a + b")
   df["a"] + df["b"]  # same as the previous expression

In certain cases :meth:`~pandas.DataFrame.eval` will be much faster than
evaluation in pure Python. For more details and examples see :ref:`the eval
documentation <enhancingperf.eval>`.

plyr
----

``plyr`` is an R library for the split-apply-combine strategy for data
analysis. The functions revolve around three data structures in R, ``a``
for ``arrays``, ``l`` for ``lists``, and ``d`` for ``data.frame``. The
table below shows how these data structures could be mapped in Python.

+------------+-------------------------------+
| R          | Python                        |
+============+===============================+
| array      | list                          |
+------------+-------------------------------+
| lists      | dictionary or list of objects |
+------------+-------------------------------+
| data.frame | dataframe                     |
+------------+-------------------------------+

|ddply|_
~~~~~~~~

An expression using a data.frame called ``df`` in R where you want to
summarize ``x`` by ``month``:

.. code-block:: r

   require(plyr)
   df <- data.frame(
     x = runif(120, 1, 168),
     y = runif(120, 7, 334),
     z = runif(120, 1.7, 20.7),
     month = rep(c(5,6,7,8),30),
     week = sample(1:4, 120, TRUE)
   )

   ddply(df, .(month, week), summarize,
         mean = round(mean(x), 2),
         sd = round(sd(x), 2))

In pandas the equivalent expression, using the
:meth:`~pandas.DataFrame.groupby` method, would be:

.. ipython:: python

   df = pd.DataFrame(
       {
           "x": np.random.uniform(1.0, 168.0, 120),
           "y": np.random.uniform(7.0, 334.0, 120),
           "z": np.random.uniform(1.7, 20.7, 120),
           "month": [5, 6, 7, 8] * 30,
           "week": np.random.randint(1, 4, 120),
       }
   )

   grouped = df.groupby(["month", "week"])
   grouped["x"].agg([np.mean, np.std])


For more details and examples see :ref:`the groupby documentation
<groupby.aggregate>`.

reshape / reshape2
------------------

|meltarray|_
~~~~~~~~~~~~~

An expression using a 3 dimensional array called ``a`` in R where you want to
melt it into a data.frame:

.. code-block:: r

   a <- array(c(1:23, NA), c(2,3,4))
   data.frame(melt(a))

In Python, since ``a`` is a list, you can simply use list comprehension.

.. ipython:: python

   a = np.array(list(range(1, 24)) + [np.NAN]).reshape(2, 3, 4)
   pd.DataFrame([tuple(list(x) + [val]) for x, val in np.ndenumerate(a)])

|meltlist|_
~~~~~~~~~~~~

An expression using a list called ``a`` in R where you want to melt it
into a data.frame:

.. code-block:: r

   a <- as.list(c(1:4, NA))
   data.frame(melt(a))

In Python, this list would be a list of tuples, so
:meth:`~pandas.DataFrame` method would convert it to a dataframe as required.

.. ipython:: python

   a = list(enumerate(list(range(1, 5)) + [np.NAN]))
   pd.DataFrame(a)

For more details and examples see :ref:`the Into to Data Structures
documentation <dsintro>`.

|meltdf|_
~~~~~~~~~~~~~~~~

An expression using a data.frame called ``cheese`` in R where you want to
reshape the data.frame:

.. code-block:: r

   cheese <- data.frame(
     first = c('John', 'Mary'),
     last = c('Doe', 'Bo'),
     height = c(5.5, 6.0),
     weight = c(130, 150)
   )
   melt(cheese, id=c("first", "last"))

In Python, the :meth:`~pandas.melt` method is the R equivalent:

.. ipython:: python

   cheese = pd.DataFrame(
       {
           "first": ["John", "Mary"],
           "last": ["Doe", "Bo"],
           "height": [5.5, 6.0],
           "weight": [130, 150],
       }
   )

   pd.melt(cheese, id_vars=["first", "last"])
   cheese.set_index(["first", "last"]).stack()  # alternative way

For more details and examples see :ref:`the reshaping documentation
<reshaping.melt>`.

|cast|_
~~~~~~~

In R ``acast`` is an expression using a data.frame called ``df`` in R to cast
into a higher dimensional array:

.. code-block:: r

   df <- data.frame(
     x = runif(12, 1, 168),
     y = runif(12, 7, 334),
     z = runif(12, 1.7, 20.7),
     month = rep(c(5,6,7),4),
     week = rep(c(1,2), 6)
   )

   mdf <- melt(df, id=c("month", "week"))
   acast(mdf, week ~ month ~ variable, mean)

In Python the best way is to make use of :meth:`~pandas.pivot_table`:

.. ipython:: python

   df = pd.DataFrame(
       {
           "x": np.random.uniform(1.0, 168.0, 12),
           "y": np.random.uniform(7.0, 334.0, 12),
           "z": np.random.uniform(1.7, 20.7, 12),
           "month": [5, 6, 7] * 4,
           "week": [1, 2] * 6,
       }
   )

   mdf = pd.melt(df, id_vars=["month", "week"])
   pd.pivot_table(
       mdf,
       values="value",
       index=["variable", "week"],
       columns=["month"],
       aggfunc=np.mean,
   )

Similarly for ``dcast`` which uses a data.frame called ``df`` in R to
aggregate information based on ``Animal`` and ``FeedType``:

.. code-block:: r

   df <- data.frame(
     Animal = c('Animal1', 'Animal2', 'Animal3', 'Animal2', 'Animal1',
                'Animal2', 'Animal3'),
     FeedType = c('A', 'B', 'A', 'A', 'B', 'B', 'A'),
     Amount = c(10, 7, 4, 2, 5, 6, 2)
   )

   dcast(df, Animal ~ FeedType, sum, fill=NaN)
   # Alternative method using base R
   with(df, tapply(Amount, list(Animal, FeedType), sum))

Python can approach this in two different ways. Firstly, similar to above
using :meth:`~pandas.pivot_table`:

.. ipython:: python

   df = pd.DataFrame(
       {
           "Animal": [
               "Animal1",
               "Animal2",
               "Animal3",
               "Animal2",
               "Animal1",
               "Animal2",
               "Animal3",
           ],
           "FeedType": ["A", "B", "A", "A", "B", "B", "A"],
           "Amount": [10, 7, 4, 2, 5, 6, 2],
       }
   )

   df.pivot_table(values="Amount", index="Animal", columns="FeedType", aggfunc="sum")

The second approach is to use the :meth:`~pandas.DataFrame.groupby` method:

.. ipython:: python

   df.groupby(["Animal", "FeedType"])["Amount"].sum()

For more details and examples see :ref:`the reshaping documentation
<reshaping.pivot>` or :ref:`the groupby documentation<groupby.split>`.

|factor|_
~~~~~~~~~

pandas has a data type for categorical data.

.. code-block:: r

   cut(c(1,2,3,4,5,6), 3)
   factor(c(1,2,3,2,2,3))

In pandas this is accomplished with ``pd.cut`` and ``astype("category")``:

.. ipython:: python

   pd.cut(pd.Series([1, 2, 3, 4, 5, 6]), 3)
   pd.Series([1, 2, 3, 2, 2, 3]).astype("category")

For more details and examples see :ref:`categorical introduction <categorical>` and the
:ref:`API documentation <api.arrays.categorical>`. There is also a documentation regarding the
:ref:`differences to R's factor <categorical.rfactor>`.


.. |c| replace:: ``c``
.. _c: https://stat.ethz.ch/R-manual/R-patched/library/base/html/c.html

.. |aggregate| replace:: ``aggregate``
.. _aggregate: https://stat.ethz.ch/R-manual/R-patched/library/stats/html/aggregate.html

.. |match| replace:: ``match`` / ``%in%``
.. _match: https://stat.ethz.ch/R-manual/R-patched/library/base/html/match.html

.. |tapply| replace:: ``tapply``
.. _tapply: https://stat.ethz.ch/R-manual/R-patched/library/base/html/tapply.html

.. |with| replace:: ``with``
.. _with: https://stat.ethz.ch/R-manual/R-patched/library/base/html/with.html

.. |subset| replace:: ``subset``
.. _subset: https://stat.ethz.ch/R-manual/R-patched/library/base/html/subset.html

.. |ddply| replace:: ``ddply``
.. _ddply: https://cran.r-project.org/web/packages/plyr/plyr.pdf#Rfn.ddply.1

.. |meltarray| replace:: ``melt.array``
.. _meltarray: https://cran.r-project.org/web/packages/reshape2/reshape2.pdf#Rfn.melt.array.1

.. |meltlist| replace:: ``melt.list``
.. meltlist: https://cran.r-project.org/web/packages/reshape2/reshape2.pdf#Rfn.melt.list.1

.. |meltdf| replace:: ``melt.data.frame``
.. meltdf: https://cran.r-project.org/web/packages/reshape2/reshape2.pdf#Rfn.melt.data.frame.1

.. |cast| replace:: ``cast``
.. cast: https://cran.r-project.org/web/packages/reshape2/reshape2.pdf#Rfn.cast.1

.. |factor| replace:: ``factor``
.. _factor: https://stat.ethz.ch/R-manual/R-devel/library/base/html/factor.html
