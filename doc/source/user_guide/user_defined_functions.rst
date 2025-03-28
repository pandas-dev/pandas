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

User-Defined Functions can be applied across various pandas methods that work with Series and DataFrames:

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
    df_applied = df.apply(add_one)
    print(df_applied)

    # This works with lambda functions too
    df_lambda = df.apply(lambda x : x + 1)
    print(df_lambda)


:meth:`DataFrame.apply` also accepts dictionaries of multiple user-defined functions:

.. ipython:: python
    # Sample DataFrame
    df = pd.DataFrame({'A': [1, 2, 3], 'B': [1, 2, 3]})
    
    # User-Defined Function
    def add_one(x):
        return x + 1

    def add_two(x):
        return x + 2
    
    # Apply function
    df_applied = df.apply({"A": add_one, "B": add_two})
    print(df_applied)

    # This works with lambda functions too
    df_lambda = df.apply({"A": lambda x : x + 1, "B": lambda x : x + 2})
    print(df_lambda)

:meth:`DataFrame.apply` works with Series objects as well:

.. ipython:: python
    # Sample Series
    s = pd.Series([1, 2, 3])
    
    # User-Defined Function
    def add_one(x):
        return x + 1
    
    # Apply function
    s_applied = s.apply(add_one)
    print(s_applied)

    # This works with lambda functions too
    s_lambda = s.apply(lambda x : x + 1)
    print(s_lambda)

:meth:`DataFrame.agg`
---------------------

The :meth:`DataFrame.agg` allows aggregation with a user-defined function along either axis (rows or columns):

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

In terms of the API, :meth:`DataFrame.agg` has similar usage to :meth:`DataFrame.apply`,
but it is primarily used for **aggregation**, applying functions that summarize or reduce data.
Typically, the result of :meth:`DataFrame.agg` reduces the dimensions of data as shown
in the above example. Conversely, :meth:`DataFrame.apply` is more general and allows for both
transformations and custom row-wise or element-wise operations.

:meth:`DataFrame.transform`
---------------------------

The :meth:`DataFrame.transform` allows transforms a Dataframe, Series or Grouped object
while preserving the original shape of the object.

.. ipython:: python 
    # Sample DataFrame  
    df = pd.DataFrame({'A': [1, 2, 3], 'B': [4, 5, 6]})  

    # User-Defined Function  
    def double(x):  
        return x * 2  

    # Apply transform  
    df_transformed = df.transform(double)  
    print(df_transformed)  

    # This works with lambda functions too  
    df_lambda = df.transform(lambda x: x * 2)  
    print(df_lambda)  

Attempting to use common aggregation functions such as `mean` or `sum` will result in
values being broadcasted to the original dimensions:

.. ipython:: python 
    # Sample DataFrame
    df = pd.DataFrame({
        'Category': ['A', 'A', 'B', 'B', 'B'],
        'Values': [10, 20, 30, 40, 50]
    })

    # Using transform with mean
    df['Mean_Transformed'] = df.groupby('Category')['Values'].transform('mean')

    # Using transform with sum
    df['Sum_Transformed'] = df.groupby('Category')['Values'].transform('sum')

    # Result broadcasted to DataFrame
    print(df)

:meth:`DataFrame.filter`
------------------------

The :meth:`DataFrame.filter` method is used to select subsets of the DataFrame’s
columns or rows and accepts user-defined functions. Specifically, these functions
return boolean values to filter columns or rows. It is useful when you want to 
extract specific columns or rows that match particular conditions.

.. ipython:: python 
    # Sample DataFrame
    df = pd.DataFrame({
        'A': [1, 2, 3],
        'B': [4, 5, 6],
        'C': [7, 8, 9],
        'D': [10, 11, 12]
    })

    # Define a function that filters out columns where the name is longer than 1 character
    df_filtered_func = df.filter(items=lambda x: len(x) > 1)
    print(df_filtered_func)

Unlike the methods discussed earlier, :meth:`DataFrame.filter` does not accept
functions that do not return boolean values, such as `mean` or `sum`.


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