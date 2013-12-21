.. currentmodule:: pandas
.. _compare_with_r:

.. ipython:: python
   :suppress:

   from pandas import *
   import numpy.random as random
   from numpy import *
   options.display.max_rows=15

Comparison with R / R libraries
*******************************

Since ``pandas`` aims to provide a lot of the data manipulation and analysis
functionality that people use `R <http://www.r-project.org/>`__ for, this page
was started to provide a more detailed look at the `R language
<http://en.wikipedia.org/wiki/R_(programming_language)>`__ and its many third
party libraries as they relate to ``pandas``. In comparisons with R and CRAN
libraries, we care about the following things:

  - **Functionality / flexibility**: what can/cannot be done with each tool
  - **Performance**: how fast are operations. Hard numbers/benchmarks are
    preferable
  - **Ease-of-use**: Is one tool easier/harder to use (you may have to be
    the judge of this, given side-by-side code comparisons)

This page is also here to offer a bit of a translation guide for users of these
R packages.

Base R
------

|subset|_
~~~~~~~~~~

.. versionadded:: 0.13

The :meth:`~pandas.DataFrame.query` method is similar to the base R ``subset``
function. In R you might want to get the rows of a ``data.frame`` where one
column's values are less than another column's values:

    .. code-block:: r

       df <- data.frame(a=rnorm(10), b=rnorm(10))
       subset(df, a <= b)
       df[df$a <= df$b,]  # note the comma

In ``pandas``, there are a few ways to perform subsetting. You can use
:meth:`~pandas.DataFrame.query` or pass an expression as if it were an
index/slice as well as standard boolean indexing:

    .. ipython:: python

       from pandas import DataFrame
       from numpy.random import randn

       df = DataFrame({'a': randn(10), 'b': randn(10)})
       df.query('a <= b')
       df[df.a <= df.b]
       df.loc[df.a <= df.b]

For more details and examples see :ref:`the query documentation
<indexing.query>`.


|with|_
~~~~~~~~

.. versionadded:: 0.13

An expression using a data.frame called ``df`` in R with the columns ``a`` and
``b`` would be evaluated using ``with`` like so:

    .. code-block:: r

       df <- data.frame(a=rnorm(10), b=rnorm(10))
       with(df, a + b)
       df$a + df$b  # same as the previous expression

In ``pandas`` the equivalent expression, using the
:meth:`~pandas.DataFrame.eval` method, would be:

    .. ipython:: python

       df = DataFrame({'a': randn(10), 'b': randn(10)})
       df.eval('a + b')
       df.a + df.b  # same as the previous expression

In certain cases :meth:`~pandas.DataFrame.eval` will be much faster than
evaluation in pure Python. For more details and examples see :ref:`the eval
documentation <enhancingperf.eval>`.

zoo
---

xts
---

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

In ``pandas`` the equivalent expression, using the
:meth:`~pandas.DataFrame.groupby` method, would be:



    .. ipython:: python

       df = DataFrame({
           'x': random.uniform(1., 168., 120),
           'y': random.uniform(7., 334., 120),
           'z': random.uniform(1.7, 20.7, 120),
           'month': [5,6,7,8]*30,
           'week': random.randint(1,4, 120)
       })

       grouped = df.groupby(['month','week'])
       print grouped['x'].agg([mean, std])


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
       a = array(range(1,24)+[NAN]).reshape(2,3,4)
       DataFrame([tuple(list(x)+[val]) for x, val in ndenumerate(a)])

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

       a = list(enumerate(range(1,5)+[NAN]))
       DataFrame(a)

For more details and examples see :ref:`the Into to Data Structures
documentation <basics.dataframe.from_items>`.

|meltdf|_
~~~~~~~~~~~~~~~~

An expression using a data.frame called ``cheese`` in R where you want to
reshape the data.frame:

    .. code-block:: r

       cheese <- data.frame(
         first = c('John, Mary'),
         last = c('Doe', 'Bo'),
         height = c(5.5, 6.0),
         weight = c(130, 150)
       )
       melt(cheese, id=c("first", "last"))

In Python, the :meth:`~pandas.melt` method is the R equivalent:

    .. ipython:: python

       cheese = DataFrame({'first' : ['John', 'Mary'],
                           'last' : ['Doe', 'Bo'],
                           'height' : [5.5, 6.0],
                           'weight' : [130, 150]})
       melt(cheese, id_vars=['first', 'last'])
       cheese.set_index(['first', 'last']).stack() # alternative way

For more details and examples see :ref:`the reshaping documentation
<reshaping.melt>`.

|cast|_
~~~~~~~

An expression using a data.frame called ``df`` in R to cast into a higher
dimensional array:

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

        df = DataFrame({
            'x': random.uniform(1., 168., 12),
            'y': random.uniform(7., 334., 12),
            'z': random.uniform(1.7, 20.7, 12),
            'month': [5,6,7]*4,
            'week': [1,2]*6
        })
        mdf = melt(df, id_vars=['month', 'week'])
        pivot_table(mdf, values='value', rows=['variable','week'],
                    cols=['month'], aggfunc=mean)

For more details and examples see :ref:`the reshaping documentation
<reshaping.pivot>`.

.. |with| replace:: ``with``
.. _with: http://finzi.psych.upenn.edu/R/library/base/html/with.html

.. |subset| replace:: ``subset``
.. _subset: http://finzi.psych.upenn.edu/R/library/base/html/subset.html

.. |ddply| replace:: ``ddply``
.. _ddply: http://www.inside-r.org/packages/cran/plyr/docs/ddply

.. |meltarray| replace:: ``melt.array``
.. _meltarray: http://www.inside-r.org/packages/cran/reshape2/docs/melt.array

.. |meltlist| replace:: ``melt.list``
.. meltlist: http://www.inside-r.org/packages/cran/reshape2/docs/melt.list

.. |meltdf| replace:: ``melt.data.frame``
.. meltdf: http://www.inside-r.org/packages/cran/reshape2/docs/melt.data.frame

.. |cast| replace:: ``cast``
.. cast: http://www.inside-r.org/packages/cran/reshape2/docs/cast

