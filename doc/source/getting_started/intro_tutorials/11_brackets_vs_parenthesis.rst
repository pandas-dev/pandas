.. _10min_tut_11_brackets_vs_parenthesis:

Parentheses vs. Square Brackets in Python and pandas
====================================================

Python and pandas use both **parentheses** ``()`` and **square brackets** ``[]``, but these have separate and important roles. Understanding their differences is essential for writing correct Python and pandas code.

Overview of Bracket Usage
-------------------------

+------------------------------------+-------------------------------+------------------------------+
| Context                            | Parentheses ``()``            | Square Brackets ``[]``       |
+====================================+===============================+==============================+
| Function and method calls          | Yes: ``df.mean()``            | No                           |
+------------------------------------+-------------------------------+------------------------------+
| Tuple creation                     | Yes: ``(a, b, c)``            | No                           |
+------------------------------------+-------------------------------+------------------------------+
| List creation and item access      | No                            | Yes: ``a[0]``, ``[1, 2]``    |
+------------------------------------+-------------------------------+------------------------------+
| Dictionary item access             | No                            | Yes: ``mydict['key']``       |
+------------------------------------+-------------------------------+------------------------------+
| pandas column selection            | No                            | Yes: ``df['A']``,            |
|                                    |                               |     ``df[['A', 'B']]``       |
+------------------------------------+-------------------------------+------------------------------+
| pandas row selection and slicing   | No                            | Yes: ``df[0:5]``,            |
|                                    |                               |     ``df.loc[2]``            |
+------------------------------------+-------------------------------+------------------------------+
| Boolean indexing                   | No                            | Yes: ``df[df['col'] > 10]``  |
+------------------------------------+-------------------------------+------------------------------+
| Grouping/logical expressions       | Yes: ``(x + y) * z``          | No                           |
+------------------------------------+-------------------------------+------------------------------+

Detailed Explanations and Examples
----------------------------------

**Parentheses ``()``**

- Used to *call* functions and methods, enclosing arguments:
  .. code-block:: python

      df.mean()
      print("Hello, world!")
      sum([1, 2, 3])

- Used to create *tuples*, which are immutable ordered collections:
  .. code-block:: python

      coordinates = (10, 20)
      empty_tuple = ()

- Used for *grouping expressions* in mathematics and logic:
  .. code-block:: python

      result = (1 + 2) * 3
      is_valid = (score > 0) and (score < 100)

- Used to spread Python statements over multiple lines:

  .. code-block:: python

      total = (
          1 +
          2 +
          3
      )

**Square Brackets ``[]``**

- Used to *define* and *access* elements of Python lists:
  .. code-block:: python

      numbers = [1, 2, 3, 4]
      first = numbers[0]
      sub = numbers[1:3]

- Used to *access* values in dictionaries by key:
  .. code-block:: python

      prices = {'apple': 40, 'banana': 10}
      apple_price = prices['apple']

- Used for all kinds of *indexing* and *selection* in pandas DataFrames and Series:

  *Single column selection* (returns Series):
  .. code-block:: python

      df['A']

  *Multiple columns selection* (returns DataFrame):
  .. code-block:: python

      df[['A', 'B']]

    Here, the inner brackets create a Python list of column labels, and the outer brackets are pandas selection syntax.

  *Row selection and slicing*:
  .. code-block:: python

      df[0:2]           # selects rows 0 and 1 by integer position
      df.loc[2]         # selects row with label/index 2
      df.iloc[2]        # selects row at integer position 2

  *Boolean indexing (row filtering)*:
  .. code-block:: python

      df[df['A'] > 5]   # returns only rows where column 'A' is greater than 5

Key Points to Remember
----------------------

- **Parentheses** are for function/method calls, tuple creation, grouping, and continuation of statements.
- **Square brackets** are for creating and accessing lists, dictionary values, slicing sequences, and—critically for pandas—selecting/subsetting columns and rows.
- In pandas, *single* square brackets select a single column as a Series (``df['A']``), while *double* square brackets select multiple columns as a DataFrame (``df[['A', 'B']]``) because the *inner brackets create a Python list* of column labels.
- Boolean indexing in pandas always uses square brackets enclosing a boolean Series: ``df[df['A'] > 5]``.

Common Pitfalls
---------------

- Attempting to use parentheses for list/tensor/column access or slicing will result in errors.
- Using single brackets with a list inside (like ``df[['A']]``) still returns a DataFrame, not a Series—bracket count and the type of object inside matters.
- Remember that method calls (like ``df.mean()`` or ``df.groupby('A')``) always require parentheses, even if there are no arguments.

Example Summary
---------------

.. code-block:: python

    import pandas as pd
    df = pd.DataFrame({'A': [1, 2, 3], 'B': [4, 5, 6]})

    # Function and method calls
    df.mean()

    # Tuple creation
    t = (1, 2, 3)

    # List creation/access
    mylist = [10, 20, 30]
    first_item = mylist[0]

    # pandas column selection
    df['A']
    df[['A', 'B']]

    # pandas boolean indexing
    df[df['B'] > 4]

    # Grouping/logical expressions
    selected = (df['A'] > 1) & (df['B'] < 6)

Getting comfortable with the distinction between parentheses and square brackets is a major milestone for every new pandas user. This understanding leads to correct code and enhanced productivity when working in both core Python and pandas.
