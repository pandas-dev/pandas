.. currentmodule:: pandas
.. _compare_with_r:

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

reshape / reshape2
------------------


.. |with| replace:: ``with``
.. _with: http://finzi.psych.upenn.edu/R/library/base/html/with.html

.. |subset| replace:: ``subset``
.. _subset: http://finzi.psych.upenn.edu/R/library/base/html/subset.html
