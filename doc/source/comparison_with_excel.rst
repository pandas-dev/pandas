.. currentmodule:: pandas 
.. _compare_with_excel:

.. ipython:: python
    :suppress:

    import pandas as pd 
    import random
    pd.options.display.max_rows=15

Comparison with Excel
*********************

Commonly used Excel functionalities
-----------------------------------

Fill Handle
~~~~~~~~~~~

Create a series of numbers following a set pattern in a certain set of cells. In
Excel this would be done by shift+drag after entering the first number or by 
entering the first two or three values and then dragging.

This can be achieved by creating a series and assigning it to the desired cells.

.. ipython:: python

    df = pd.DataFrame({'AAA': [1] * 8, 'BBB': list(range(0, 8))}); df

    series = list(range(1, 5)); series

    df.iloc[2:(5+1)].AAA = series

    df

Filters
~~~~~~~

Filters can be achieved by using slicing. 

The examples filter by 0 on column AAA, and also show how to filter by multiple
values.

.. ipython:: python

   df[df.AAA == 0]

   df[(df.AAA == 0) | (df.AAA == 2)]


Drop Duplicates
~~~~~~~~~~~~~~~

Another commonly used function is Drop Duplicates. This is directly supported in 
pandas.

.. ipython:: python

    df = pd.DataFrame({"class": ['A', 'A', 'A', 'B', 'C', 'D'], "student_count": [42, 35, 42, 50, 47, 45], "all_pass": ["Yes", "Yes", "Yes", "No", "No", "Yes"]})

    df.drop_duplicates()

    df.drop_duplicates(["class", "student_count"])


Pivot Table
~~~~~~~~~~~

This can be achieved by using ``pandas.pivot_table`` for examples and reference, 
please see `pandas.pivot_table <http://pandas.pydata.org/pandas-docs/stable/generated/pandas.pivot_table.html>`__ 


Formulae
~~~~~~~~

Let's create a new column "girls_count" and try to compute the number of boys in 
each class. 

.. ipython:: python 

    df["girls_count"]  = [21, 12, 21, 31, 23, 17]; df

    def get_count(row):
        return row["student_count"] - row["girls_count"]

    df["boys_count"] = df.apply(get_count, axis = 1); df


VLOOKUP
~~~~~~~

.. ipython:: python

    df1 = pd.DataFrame({"keys": [1, 2, 3, 4, 5, 6, 7], "first_names": ["harry", "ron",
    "hermione", "rubius", "albus", "severus", "luna"]}); df1

    random_names = pd.DataFrame({"surnames": ["hadrid", "malfoy", "lovegood",
    "dumbledore", "grindelwald", "granger", "weasly", "riddle", "longbottom",
    "snape"], "keys": [ random.randint(1,7) for x in range(0,10) ]})

    random_names

    random_names.merge(df1, on="keys", how='left')

Adding a row
~~~~~~~~~~~~

To appended a row, we can just assign values to an index using ``iloc``. 

NOTE: If the index already exists, the values in that index will be over written.

.. ipython:: python

    df1.iloc[7] = [8, "tonks"]; df1


Search and Replace
~~~~~~~~~~~~~~~~~~

The ``replace`` method that comes associated with the ``DataFrame`` object can perform 
this function. Please see `pandas.DataFrame.replace <https://pandas.pydata.org/pandas-docs/stable/generated/pandas.DataFrame.replace.html>`__ for examples.
