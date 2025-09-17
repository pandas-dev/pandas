Parentheses vs. Brackets in pandas
==================================

In pandas, beginners often get confused between **parentheses ``()``** and **square brackets ``[]``**.
Understanding the difference is essential for using pandas effectively.

**Parentheses ``()``** are used to **call functions or methods**:

.. code-block:: python

    import pandas as pd
    df = pd.DataFrame({'A': [1, 2, 3], 'B': [4, 5, 6]})

    # Call the head() method to see first 5 rows
    df.head()

**Square brackets ``[]``** are used to **access data or select columns**:

.. code-block:: python

    # Select a single column as a Series
    df['A']

    # Select multiple columns as a DataFrame
    df[['A', 'B']]

**Key points:**

- `()` always executes something (a function/method).
- `[]` always retrieves data (like indexing or slicing).
- Mixing them up is a common source of errors for new pandas users.

**Additional examples:**

.. code-block:: python

    # Using brackets to filter rows
    df[df['A'] > 1]

    # Using parentheses to chain method calls
    df[['A', 'B']].head()
