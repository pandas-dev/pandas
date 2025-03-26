.. _user_defined_functions:

{{ header }}

**************************************
Introduction to User-Defined Functions
**************************************

In pandas, User-Defined Functions (UDFs) provide a way to extend the library’s
functionality by allowing users to apply custom computations to their data. While
pandas comes with a set of built-in functions for data manipulation, UDFs offer
flexibility when built-in methods are not sufficient. These functions can be 
applied at different levels: element-wise, row-wise, column-wise, or group-wise,
depending on the method used.

.. .. note:: 
    
..     User-Defined Functions will be abbreviated to UDFs throughout this guide.

Why Use User-Defined Functions?
-------------------------------

Pandas is designed for high-performance data processing, but sometimes your specific
needs go beyond standard aggregation, transformation, or filtering. UDFs allow you to:

* **Customize Computations**: Implement logic tailored to your dataset, such as complex 
  transformations, domain-specific calculations, or conditional modifications.
* **Improve Code Readability**: Encapsulate logic into functions rather than writing long,
  complex expressions.
* **Handle Complex Grouped Operations**: Perform operations on grouped data that standard
  methods do not support.
* **Extend pandas' Functionality**: Apply external libraries or advanced calculations that 
  are not natively available.


What functions support User-Defined Functions
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

UDFs can be applied across various pandas methods that work with Series and DataFrames:

* :meth:`DataFrame.apply` - A flexible method that allows applying a function to Series,
  DataFrames, or groups of data.
* :meth:`DataFrame.agg` (Aggregate) - Used for summarizing data, supporting multiple
  aggregation functions.
* :meth:`DataFrame.transform` - Applies a function to groups while preserving the shape of
  the original data.
* :meth:`DataFrame.filter` - Filters groups based on a function returning a Boolean condition.
* :meth:`DataFrame.map` - Applies an element-wise function to a Series, useful for
  transforming individual values.
* :meth:`DataFrame.pipe` - Allows chaining custom functions to process entire DataFrames or
  Series in a clean, readable manner.

Each of these methods can be used with both Series and DataFrame objects, providing versatile
ways to apply user-defined functions across different pandas data structures.


:meth:`DataFrame.apply`
-----------------------

The :meth:`DataFrame.apply` allows applying a user-defined functions along either axis (rows or columns):

.. ipython:: python

    import pandas as pd
    
    # Sample DataFrame
    df = pd.DataFrame({'A': [1, 2, 3], 'B': [4, 5, 6]})
    
    # User-Defined Function
    def add_one(x):
        return x + 1
    
    # Apply function
    df_transformed = df.apply(add_one)
    print(df_transformed)

    # This works with lambda functions too
    df_lambda = df.apply(lambda x : x + 1)
    print(df_lambda)


:meth:`DataFrame.apply` also accepts dictionaries of multiple user-defined functions:

.. ipython:: python

    import pandas as pd
    
    # Sample DataFrame
    df = pd.DataFrame({'A': [1, 2, 3], 'B': [1, 2, 3]})
    
    # User-Defined Function
    def add_one(x):
        return x + 1

    def add_two(x):
        return x + 2
    
    # Apply function
    df_transformed = df.apply({"A": add_one, "B": add_two})
    print(df_transformed)

    # This works with lambda functions too
    df_lambda = df.apply({"A": lambda x : x + 1, "B": lambda x : x + 2})
    print(df_lambda)

:meth:`DataFrame.apply` works with Series objects as well:

.. ipython:: python

    import pandas as pd
    
    # Sample Series
    s = pd.Series([1, 2, 3])
    
    # User-Defined Function
    def add_one(x):
        return x + 1
    
    # Apply function
    s_transformed = s.apply(add_one)
    print(df_transformed)

    # This works with lambda functions too
    s_lambda = s.apply(lambda x : x + 1)
    print(s_lambda)

:meth:`DataFrame.agg`
---------------------

When working with grouped data, user-defined functions can be used within :meth:`DataFrame.agg`:

.. ipython:: python

    # Sample DataFrame
    df = pd.DataFrame({
        'Category': ['A', 'A', 'B', 'B'],
        'Values': [10, 20, 30, 40]
    })
    
    # Define a function for group operations
    def group_mean(group):
        return group.mean()
    
    # Apply UDF to each group
    grouped_result = df.groupby('Category')['Values'].agg(group_mean)
    print(grouped_result)

Performance Considerations
--------------------------

While UDFs provide flexibility, their use is currently discouraged as they can introduce performance issues, especially when
written in pure Python. To improve efficiency:

* Use **vectorized operations** (`NumPy` or `pandas` built-ins) when possible.
* Leverage **Cython or Numba** to speed up computations.
* Consider using **pandas' built-in methods** instead of UDFs for common operations.

.. note::
    If performance is critical, explore **pandas' vectorized functions** before resorting
    to UDFs.