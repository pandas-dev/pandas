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
and change the data differently, depending on the method used.

Why Use User-Defined Functions?
-------------------------------

Pandas is designed for high-performance data processing, but sometimes your specific
needs go beyond standard aggregation, transformation, or filtering. User-defined functions allow you to:

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

User-Defined Functions can be applied across various pandas methods:

* :meth:`DataFrame.apply` - A flexible method that allows applying a function to Series,
  DataFrames, or groups of data.
* :meth:`DataFrame.agg` (Aggregate) - Used for summarizing data, supporting multiple
  aggregation functions.
* :meth:`DataFrame.transform` - Applies a function to groups while preserving the shape of
  the original data.
* :meth:`DataFrame.filter` - Filters groups based on a list of Boolean conditions.
* :meth:`DataFrame.map` - Applies an element-wise function to a Series, useful for
  transforming individual values.
* :meth:`DataFrame.pipe` - Allows chaining custom functions to process entire DataFrames or
  Series in a clean, readable manner.

All of these pandas methods can be used with both Series and DataFrame objects, providing versatile
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
columns or row. It is useful when you want to extract specific columns or rows that
match particular conditions.

.. note::
    :meth:`DataFrame.filter` does not accept user-defined functions, but can accept
    list comprehensions that have user-defined functions applied to them.

.. ipython:: python 

    # Sample DataFrame
    df = pd.DataFrame({
        'AA': [1, 2, 3],
        'BB': [4, 5, 6],
        'C': [7, 8, 9],
        'D': [10, 11, 12]
    })

    def is_long_name(column_name):
        return len(column_name) > 1

    # Define a function that filters out columns where the name is longer than 1 character
    df_filtered = df[[col for col in df.columns if is_long_name(col)]]
    print(df_filtered)

:meth:`DataFrame.map`
---------------------

The :meth:`DataFrame.map` method is used to apply a function element-wise to a pandas Series
or Dataframe. It is particularly useful for substituting values or transforming data.

.. ipython:: python

    # Sample DataFrame
    df = pd.DataFrame({ 'A': ['cat', 'dog', 'bird'], 'B': ['pig', 'cow', 'lamb'] })

    # Using map with a user-defined function
    def animal_to_length(animal):
        return len(animal)

    df_mapped = df.map(animal_to_length)
    print(df_mapped)

    # This works with lambda functions too
    df_lambda = df.map(lambda x: x.upper())
    print(df_lambda)

:meth:`DataFrame.pipe`
----------------------

The :meth:`DataFrame.pipe` method allows you to apply a function or a series of functions to a
DataFrame in a clean and readable way. This is especially useful for building data processing pipelines.

.. ipython:: python

    # Sample DataFrame
    df = pd.DataFrame({ 'A': [1, 2, 3], 'B': [4, 5, 6] })

    # User-defined functions for transformation
    def add_one(df):
        return df + 1

    def square(df):
        return df ** 2

    # Applying functions using pipe
    df_piped = df.pipe(add_one).pipe(square)
    print(df_piped)

The advantage of using :meth:`DataFrame.pipe` is that it allows you to chain together functions
without nested calls, promoting a cleaner and more readable code style.


Performance Considerations
--------------------------

While user-defined functions provide flexibility, their use is currently discouraged as they can introduce
performance issues, especially when written in pure Python. To improve efficiency,
consider using built-in `NumPy` or `pandas` functions instead of user-defined functions
for common operations.

.. note::
    If performance is critical, explore **vectorizated operations** before resorting
    to user-defined functions.

Vectorized Operations
~~~~~~~~~~~~~~~~~~~~~

Below is an example of vectorized operations in pandas:

.. code-block:: text

    # User-defined function
    def calc_ratio(row):
        return 100 * (row["one"] / row["two"])

    df["new_col2"] = df.apply(calc_ratio, axis=1)

    # Vectorized Operation
    df["new_col"] = 100 * (df["one"] / df["two"])

Measuring how long each operation takes:

.. code-block:: text

    Vectorized:             0.0043 secs
    User-defined function:  5.6435 secs

Vectorized operations in pandas are significantly faster than using :meth:`DataFrame.apply`
with user-defined functions because they leverage highly optimized C functions
via NumPy to process entire arrays at once. This approach avoids the overhead of looping
through rows in Python and making separate function calls for each row, which is slow and
inefficient. Additionally, NumPy arrays benefit from memory efficiency and CPU-level
optimizations, making vectorized operations the preferred choice whenever possible.