.. _compare_with_sql:

{{ header }}

Comparison with SQL
********************
Since many potential pandas users have some familiarity with
`SQL <https://en.wikipedia.org/wiki/SQL>`_, this page is meant to provide some examples of how
various SQL operations would be performed using pandas.

.. include:: includes/introduction.rst

Most of the examples will utilize the ``tips`` dataset found within pandas tests.  We'll read
the data into a DataFrame called ``tips`` and assume we have a database table of the same name and
structure.

.. ipython:: python

    url = (
        "https://raw.githubusercontent.com/pandas-dev"
        "/pandas/main/pandas/tests/io/data/csv/tips.csv"
    )
    tips = pd.read_csv(url)
    tips


Copies vs. in place operations
------------------------------

.. include:: includes/copies.rst


SELECT
------
In SQL, selection is done using a comma-separated list of columns you'd like to select (or a ``*``
to select all columns):

.. code-block:: sql

    SELECT total_bill, tip, smoker, time
    FROM tips;

With pandas, column selection is done by passing a list of column names to your DataFrame:

.. ipython:: python

    tips[["total_bill", "tip", "smoker", "time"]]

Calling the DataFrame without the list of column names would display all columns (akin to SQL's
``*``).

In SQL, you can add a calculated column:

.. code-block:: sql

    SELECT *, tip/total_bill as tip_rate
    FROM tips;

With pandas, you can use the :meth:`DataFrame.assign` method of a DataFrame to append a new column:

.. ipython:: python

    tips.assign(tip_rate=tips["tip"] / tips["total_bill"])

WHERE
-----
Filtering in SQL is done via a WHERE clause.

.. code-block:: sql

    SELECT *
    FROM tips
    WHERE time = 'Dinner';

.. include:: includes/filtering.rst

Just like SQL's ``OR`` and ``AND``, multiple conditions can be passed to a DataFrame using ``|``
(``OR``) and ``&`` (``AND``).

Tips of more than $5 at Dinner meals:

.. code-block:: sql

    SELECT *
    FROM tips
    WHERE time = 'Dinner' AND tip > 5.00;

.. ipython:: python

    tips[(tips["time"] == "Dinner") & (tips["tip"] > 5.00)]

Tips by parties of at least 5 diners OR bill total was more than $45:

.. code-block:: sql

    SELECT *
    FROM tips
    WHERE size >= 5 OR total_bill > 45;

.. ipython:: python

    tips[(tips["size"] >= 5) | (tips["total_bill"] > 45)]

NULL checking is done using the :meth:`~pandas.Series.notna` and :meth:`~pandas.Series.isna`
methods.

.. ipython:: python

    frame = pd.DataFrame(
        {"col1": ["A", "B", np.nan, "C", "D"], "col2": ["F", np.nan, "G", "H", "I"]}
    )
    frame

Assume we have a table of the same structure as our DataFrame above. We can see only the records
where ``col2`` IS NULL with the following query:

.. code-block:: sql

    SELECT *
    FROM frame
    WHERE col2 IS NULL;

.. ipython:: python

    frame[frame["col2"].isna()]

Getting items where ``col1`` IS NOT NULL can be done with :meth:`~pandas.Series.notna`.

.. code-block:: sql

    SELECT *
    FROM frame
    WHERE col1 IS NOT NULL;

.. ipython:: python

    frame[frame["col1"].notna()]


GROUP BY
--------
In pandas, SQL's ``GROUP BY`` operations are performed using the similarly named
:meth:`~pandas.DataFrame.groupby` method. :meth:`~pandas.DataFrame.groupby` typically refers to a
process where we'd like to split a dataset into groups, apply some function (typically aggregation)
, and then combine the groups together.

A common SQL operation would be getting the count of records in each group throughout a dataset.
For instance, a query getting us the number of tips left by sex:

.. code-block:: sql

    SELECT sex, count(*)
    FROM tips
    GROUP BY sex;
    /*
    Female     87
    Male      157
    */


The pandas equivalent would be:

.. ipython:: python

    tips.groupby("sex").size()

Notice that in the pandas code we used :meth:`.DataFrameGroupBy.size` and not
:meth:`.DataFrameGroupBy.count`. This is because
:meth:`.DataFrameGroupBy.count` applies the function to each column, returning
the number of ``NOT NULL`` records within each.

.. ipython:: python

    tips.groupby("sex").count()

Alternatively, we could have applied the :meth:`.DataFrameGroupBy.count` method
to an individual column:

.. ipython:: python

    tips.groupby("sex")["total_bill"].count()

Multiple functions can also be applied at once. For instance, say we'd like to see how tip amount
differs by day of the week - :meth:`.DataFrameGroupBy.agg` allows you to pass a dictionary
to your grouped DataFrame, indicating which functions to apply to specific columns.

.. code-block:: sql

    SELECT day, AVG(tip), COUNT(*)
    FROM tips
    GROUP BY day;
    /*
    Fri   2.734737   19
    Sat   2.993103   87
    Sun   3.255132   76
    Thu  2.771452   62
    */

.. ipython:: python

    tips.groupby("day").agg({"tip": "mean", "day": "size"})

Grouping by more than one column is done by passing a list of columns to the
:meth:`~pandas.DataFrame.groupby` method.

.. code-block:: sql

    SELECT smoker, day, COUNT(*), AVG(tip)
    FROM tips
    GROUP BY smoker, day;
    /*
    smoker day
    No     Fri      4  2.812500
           Sat     45  3.102889
           Sun     57  3.167895
           Thu    45  2.673778
    Yes    Fri     15  2.714000
           Sat     42  2.875476
           Sun     19  3.516842
           Thu    17  3.030000
    */

.. ipython:: python

    tips.groupby(["smoker", "day"]).agg({"tip": ["size", "mean"]})

.. _compare_with_sql.join:

JOIN
----
``JOIN``\s can be performed with :meth:`~pandas.DataFrame.join` or :meth:`~pandas.merge`. By
default, :meth:`~pandas.DataFrame.join` will join the DataFrames on their indices. Each method has
parameters allowing you to specify the type of join to perform (``LEFT``, ``RIGHT``, ``INNER``,
``FULL``) or the columns to join on (column names or indices).

.. warning::

    If both key columns contain rows where the key is a null value, those
    rows will be matched against each other. This is different from usual SQL
    join behaviour and can lead to unexpected results.

.. ipython:: python

    df1 = pd.DataFrame({"key": ["A", "B", "C", "D"], "value": np.random.randn(4)})
    df2 = pd.DataFrame({"key": ["B", "D", "D", "E"], "value": np.random.randn(4)})

Assume we have two database tables of the same name and structure as our DataFrames.

Now let's go over the various types of ``JOIN``\s.

INNER JOIN
~~~~~~~~~~
.. code-block:: sql

    SELECT *
    FROM df1
    INNER JOIN df2
      ON df1.key = df2.key;

.. ipython:: python

    # merge performs an INNER JOIN by default
    pd.merge(df1, df2, on="key")

:meth:`~pandas.merge` also offers parameters for cases when you'd like to join one DataFrame's
column with another DataFrame's index.

.. ipython:: python

    indexed_df2 = df2.set_index("key")
    pd.merge(df1, indexed_df2, left_on="key", right_index=True)

LEFT OUTER JOIN
~~~~~~~~~~~~~~~

Show all records from ``df1``.

.. code-block:: sql

    SELECT *
    FROM df1
    LEFT OUTER JOIN df2
      ON df1.key = df2.key;

.. ipython:: python

    pd.merge(df1, df2, on="key", how="left")

RIGHT JOIN
~~~~~~~~~~

Show all records from ``df2``.

.. code-block:: sql

    SELECT *
    FROM df1
    RIGHT OUTER JOIN df2
      ON df1.key = df2.key;

.. ipython:: python

    pd.merge(df1, df2, on="key", how="right")

FULL JOIN
~~~~~~~~~
pandas also allows for ``FULL JOIN``\s, which display both sides of the dataset, whether or not the
joined columns find a match. As of writing, ``FULL JOIN``\s are not supported in all RDBMS (MySQL).

Show all records from both tables.

.. code-block:: sql

    SELECT *
    FROM df1
    FULL OUTER JOIN df2
      ON df1.key = df2.key;

.. ipython:: python

    pd.merge(df1, df2, on="key", how="outer")


UNION
-----

``UNION ALL`` can be performed using :meth:`~pandas.concat`.

.. ipython:: python

    df1 = pd.DataFrame(
        {"city": ["Chicago", "San Francisco", "New York City"], "rank": range(1, 4)}
    )
    df2 = pd.DataFrame(
        {"city": ["Chicago", "Boston", "Los Angeles"], "rank": [1, 4, 5]}
    )

.. code-block:: sql

    SELECT city, rank
    FROM df1
    UNION ALL
    SELECT city, rank
    FROM df2;
    /*
             city  rank
          Chicago     1
    San Francisco     2
    New York City     3
          Chicago     1
           Boston     4
      Los Angeles     5
    */

.. ipython:: python

    pd.concat([df1, df2])

SQL's ``UNION`` is similar to ``UNION ALL``, however ``UNION`` will remove duplicate rows.

.. code-block:: sql

    SELECT city, rank
    FROM df1
    UNION
    SELECT city, rank
    FROM df2;
    -- notice that there is only one Chicago record this time
    /*
             city  rank
          Chicago     1
    San Francisco     2
    New York City     3
           Boston     4
      Los Angeles     5
    */

In pandas, you can use :meth:`~pandas.concat` in conjunction with
:meth:`~pandas.DataFrame.drop_duplicates`.

.. ipython:: python

    pd.concat([df1, df2]).drop_duplicates()


LIMIT
-----

.. code-block:: sql

    SELECT * FROM tips
    LIMIT 10;

.. ipython:: python

    tips.head(10)


pandas equivalents for some SQL analytic and aggregate functions
----------------------------------------------------------------

Top n rows with offset
~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: sql

    -- MySQL
    SELECT * FROM tips
    ORDER BY tip DESC
    LIMIT 10 OFFSET 5;

.. ipython:: python

    tips.nlargest(10 + 5, columns="tip").tail(10)

Top n rows per group
~~~~~~~~~~~~~~~~~~~~

.. code-block:: sql

    -- Oracle's ROW_NUMBER() analytic function
    SELECT * FROM (
      SELECT
        t.*,
        ROW_NUMBER() OVER(PARTITION BY day ORDER BY total_bill DESC) AS rn
      FROM tips t
    )
    WHERE rn < 3
    ORDER BY day, rn;


.. ipython:: python

    (
        tips.assign(
            rn=tips.sort_values(["total_bill"], ascending=False)
            .groupby(["day"])
            .cumcount()
            + 1
        )
        .query("rn < 3")
        .sort_values(["day", "rn"])
    )

the same using ``rank(method='first')`` function

.. ipython:: python

    (
        tips.assign(
            rnk=tips.groupby(["day"])["total_bill"].rank(
                method="first", ascending=False
            )
        )
        .query("rnk < 3")
        .sort_values(["day", "rnk"])
    )

.. code-block:: sql

    -- Oracle's RANK() analytic function
    SELECT * FROM (
      SELECT
        t.*,
        RANK() OVER(PARTITION BY sex ORDER BY tip) AS rnk
      FROM tips t
      WHERE tip < 2
    )
    WHERE rnk < 3
    ORDER BY sex, rnk;

Let's find tips with (rank < 3) per gender group for (tips < 2).
Notice that when using ``rank(method='min')`` function
``rnk_min`` remains the same for the same ``tip``
(as Oracle's ``RANK()`` function)

.. ipython:: python

    (
        tips[tips["tip"] < 2]
        .assign(rnk_min=tips.groupby(["sex"])["tip"].rank(method="min"))
        .query("rnk_min < 3")
        .sort_values(["sex", "rnk_min"])
    )


UPDATE
------

.. code-block:: sql

    UPDATE tips
    SET tip = tip*2
    WHERE tip < 2;

.. ipython:: python

    tips.loc[tips["tip"] < 2, "tip"] *= 2

DELETE
------

.. code-block:: sql

    DELETE FROM tips
    WHERE tip > 9;

In pandas we select the rows that should remain instead of deleting the rows that should be removed:

.. ipython:: python

    tips = tips.loc[tips["tip"] <= 9]
