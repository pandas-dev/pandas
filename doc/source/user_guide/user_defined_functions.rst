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

Why Not To Use User-Defined Functions
-----------------------------------------

While UDFs provide flexibility, they come with significant drawbacks, primarily
related to performance. Unlike vectorized pandas operations, UDFs are slower because pandas lacks
insight into what they are computing, making it difficult to apply efficient handling or optimization
techniques. As a result, pandas resorts to less efficient processing methods that significantly
slow down computations. Additionally, relying on UDFs often sacrifices the benefits
of pandas’ built-in, optimized methods, limiting compatibility and overall performance.

.. note::
    In general, most tasks can and should be accomplished using pandas’ built-in methods or vectorized operations.

Despite their drawbacks, UDFs can be helpful when:

* **Custom Computations Are Needed**: Implementing complex logic or domain-specific calculations that pandas'
  built-in methods cannot handle.
* **Extending pandas' Functionality**: Applying external libraries or specialized algorithms unavailable in pandas.
* **Handling Complex Grouped Operations**: Performing operations on grouped data that standard methods do not support.

Methods that support User-Defined Functions
-------------------------------------------

User-Defined Functions can be applied across various pandas methods:

* :meth:`~DataFrame.apply` - A flexible method that allows applying a function to Series,
  DataFrames, or groups of data.
* :meth:`~DataFrame.agg` (Aggregate) - Used for summarizing data, supporting multiple
  aggregation functions.
* :meth:`~DataFrame.transform` - Applies a function to groups while preserving the shape of
  the original data.
* :meth:`~DataFrame.filter` - Filters groups based on a list of Boolean conditions.
* :meth:`~DataFrame.map` - Applies an element-wise function to a Series, useful for
  transforming individual values.
* :meth:`~DataFrame.pipe` - Allows chaining custom functions to process entire DataFrames or
  Series in a clean, readable manner.

All of these pandas methods can be used with both Series and DataFrame objects, providing versatile
ways to apply UDFs across different pandas data structures.

.. note::
    Some of these methods are can also be applied to Groupby Objects. Refer to :ref:`groupby`.


Choosing the Right Method
-------------------------
When applying UDFs in pandas, it is essential to select the appropriate method based
on your specific task. Each method has its strengths and is designed for different use
cases. Understanding the purpose and behavior of each method will help you make informed
decisions, ensuring more efficient and maintainable code.

Below is a table overview of all methods that accept UDFs:

+------------------+--------------------------------------+---------------------------+--------------------+---------------------------+------------------------------------------+
| Method           | Purpose                              | Supports UDFs             | Keeps Shape        | Performance               | Recommended Use Case                     |
+==================+======================================+===========================+====================+===========================+==========================================+
| :meth:`apply`    | General-purpose function             | Yes                       | Yes (when axis=1)  | Slow                      | Custom row-wise or column-wise operations|
+------------------+--------------------------------------+---------------------------+--------------------+---------------------------+------------------------------------------+
| :meth:`agg`      | Aggregation                          | Yes                       | No                 | Fast (if using built-ins) | Custom aggregation logic                 |
+------------------+--------------------------------------+---------------------------+--------------------+---------------------------+------------------------------------------+
| :meth:`transform`| Transform without reducing dimensions| Yes                       | Yes                | Fast (if vectorized)      | Broadcast element-wise transformations   |
+------------------+--------------------------------------+---------------------------+--------------------+---------------------------+------------------------------------------+
| :meth:`map`      | Element-wise mapping                 | Yes                       | Yes                | Moderate                  | Simple element-wise transformations      |
+------------------+--------------------------------------+---------------------------+--------------------+---------------------------+------------------------------------------+
| :meth:`pipe`     | Functional chaining                  | Yes                       | Yes                | Depends on function       | Building clean pipelines                 |
+------------------+--------------------------------------+---------------------------+--------------------+---------------------------+------------------------------------------+
| :meth:`filter`   | Row/Column selection                 | Not directly              | Yes                | Fast                      | Subsetting based on conditions           |
+------------------+--------------------------------------+---------------------------+--------------------+---------------------------+------------------------------------------+

:meth:`DataFrame.apply`
~~~~~~~~~~~~~~~~~~~~~~~

The :meth:`DataFrame.apply` allows you to apply UDFs along either rows or columns. While flexible,
it is slower than vectorized operations and should be used only when you need operations
that cannot be achieved with built-in pandas functions.

When to use: :meth:`DataFrame.apply` is suitable when no alternative vectorized method is available, but consider
optimizing performance with vectorized operations wherever possible.

Examples of usage can be found :ref:`here <api.dataframe.apply>`.

:meth:`DataFrame.agg`
~~~~~~~~~~~~~~~~~~~~~

If you need to aggregate data, :meth:`DataFrame.agg` is a better choice than apply because it is
specifically designed for aggregation operations.

When to use: Use :meth:`DataFrame.agg` for performing aggregations like sum, mean, or custom aggregation
functions across groups.

Examples of usage can be found :ref:`here <api.dataframe.agg>`.

:meth:`DataFrame.transform`
~~~~~~~~~~~~~~~~~~~~~~~~~~~

The transform method is ideal for performing element-wise transformations while preserving the shape of the original DataFrame.
It’s generally faster than apply because it can take advantage of pandas' internal optimizations.

When to use: When you need to perform element-wise transformations that retain the original structure of the DataFrame.

Documentation can be found :ref:`here <api.dataframe.transform>`.

Attempting to use common aggregation functions such as ``mean`` or ``sum`` will result in
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
~~~~~~~~~~~~~~~~~~~~~~~~

The :meth:`DataFrame.filter` method is used to select subsets of the DataFrame’s
columns or row. It is useful when you want to extract specific columns or rows that
match particular conditions.

When to use: Use :meth:`DataFrame.filter` when you want to use a UDF to create a subset of a DataFrame or Series

.. note::
    :meth:`DataFrame.filter` does not accept UDFs, but can accept
    list comprehensions that have UDFs applied to them.

.. ipython:: python

    # Sample DataFrame
    df = pd.DataFrame({
        'AA': [1, 2, 3],
        'BB': [4, 5, 6],
        'C': [7, 8, 9],
        'D': [10, 11, 12]
    })

    # Define a function that filters out columns where the name is longer than 1 character
    def is_long_name(column_name):
        return len(column_name) > 1

    df_filtered = df[[col for col in df.columns if is_long_name(col)]]
    print(df_filtered)

Since filter does not directly accept a UDF, you have to apply the UDF indirectly,
such as by using list comprehensions.

:meth:`DataFrame.map`
~~~~~~~~~~~~~~~~~~~~~

:meth:`DataFrame.map` is used specifically to apply element-wise UDFs and is better
for this purpose compared to :meth:`DataFrame.apply` because of its better performance.

When to use: Use map for applying element-wise UDFs to DataFrames or Series.

Documentation can be found :ref:`here <api.dataframe.map>`.

:meth:`DataFrame.pipe`
~~~~~~~~~~~~~~~~~~~~~~

The pipe method is useful for chaining operations together into a clean and readable pipeline.
It is a helpful tool for organizing complex data processing workflows.

When to use: Use pipe when you need to create a pipeline of transformations and want to keep the code readable and maintainable.

Documentation can be found :ref:`here <api.dataframe.pipe>`.


Best Practices
--------------

While UDFs provide flexibility, their use is currently discouraged as they can introduce
performance issues, especially when written in pure Python. To improve efficiency,
consider using built-in ``NumPy`` or ``pandas`` functions instead of UDFs
for common operations.

.. note::
    If performance is critical, explore **vectorizated operations** before resorting
    to UDFs.

Vectorized Operations
~~~~~~~~~~~~~~~~~~~~~

Below is a comparison of using UDFs versus using Vectorized Operations:

.. code-block:: python

    # User-defined function
    def calc_ratio(row):
        return 100 * (row["one"] / row["two"])

    df["new_col"] = df.apply(calc_ratio, axis=1)

    # Vectorized Operation
    df["new_col2"] = 100 * (df["one"] / df["two"])

Measuring how long each operation takes:

.. code-block:: text

    User-defined function:  5.6435 secs
    Vectorized:             0.0043 secs

Vectorized operations in pandas are significantly faster than using :meth:`DataFrame.apply`
with UDFs because they leverage highly optimized C functions
via NumPy to process entire arrays at once. This approach avoids the overhead of looping
through rows in Python and making separate function calls for each row, which is slow and
inefficient. Additionally, NumPy arrays benefit from memory efficiency and CPU-level
optimizations, making vectorized operations the preferred choice whenever possible.
