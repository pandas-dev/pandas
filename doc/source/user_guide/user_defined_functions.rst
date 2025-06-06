.. _udf:

{{ header }}

*****************************
User-Defined Functions (UDFs)
*****************************

In pandas, User-Defined Functions (UDFs) provide a way to extend the library’s
functionality by allowing users to apply custom computations to their data. While
pandas comes with a set of built-in functions for data manipulation, UDFs offer
flexibility when built-in methods are not sufficient. These functions can be
applied at different levels: element-wise, row-wise, column-wise, or group-wise,
and behave differently, depending on the method used.

Here’s a simple example to illustrate a UDF applied to a Series:

.. ipython:: python

    s = pd.Series([1, 2, 3])

    # Simple UDF that adds 1 to a value
    def add_one(x):
        return x + 1

    # Apply the function element-wise using .map
    s.map(add_one)

Why Not To Use User-Defined Functions
-------------------------------------

While UDFs provide flexibility, they come with significant drawbacks, primarily
related to performance and behavior. When using UDFs, pandas must perform inference
on the result, and that inference could be incorrect. Furthermore, unlike vectorized operations,
UDFs are slower because pandas can't optimize their computations, leading to
inefficient processing.

.. note::
    In general, most tasks can and should be accomplished using pandas’ built-in methods or vectorized operations.

Despite their drawbacks, UDFs can be helpful when:

* **Custom Computations Are Needed**: Implementing complex logic or domain-specific calculations that pandas'
  built-in methods cannot handle.
* **Extending pandas' Functionality**: Applying external libraries or specialized algorithms unavailable in pandas.
* **Handling Complex Grouped Operations**: Performing operations on grouped data that standard methods do not support.

For example:

.. code-block:: python

    from sklearn.linear_model import LinearRegression

    # Sample data
    df = pd.DataFrame({
        'group': ['A', 'A', 'A', 'B', 'B', 'B'],
        'x': [1, 2, 3, 1, 2, 3],
        'y': [2, 4, 6, 1, 2, 1.5]
    })

    # Function to fit a model to each group
    def fit_model(group):
        model = LinearRegression()
        model.fit(group[['x']], group['y'])
        group['y_pred'] = model.predict(group[['x']])
        return group

    result = df.groupby('group').apply(fit_model)


Methods that support User-Defined Functions
-------------------------------------------

User-Defined Functions can be applied across various pandas methods:

+-------------------------------+------------------------+--------------------------+----------------------------------------------------------------------------------------------------------------------------------------------+
| Method                        | Function Input         | Function Output          | Description                                                                                                                                  |
+===============================+========================+==========================+==============================================================================================================================================+
| :ref:`udf.map`                | Scalar                 | Scalar                   | Apply a function to each element                                                                                                             |
+-------------------------------+------------------------+--------------------------+----------------------------------------------------------------------------------------------------------------------------------------------+
| :ref:`udf.apply` (axis=0)     | Column (Series)        | Column (Series)          | Apply a function to each column                                                                                                              |
+-------------------------------+------------------------+--------------------------+----------------------------------------------------------------------------------------------------------------------------------------------+
| :ref:`udf.apply` (axis=1)     | Row (Series)           | Row (Series)             | Apply a function to each row                                                                                                                 |
+-------------------------------+------------------------+--------------------------+----------------------------------------------------------------------------------------------------------------------------------------------+
| :ref:`udf.pipe`               | Series or DataFrame    | Series or DataFrame      | Chain functions together to apply to Series or Dataframe                                                                                     |
+-------------------------------+------------------------+--------------------------+----------------------------------------------------------------------------------------------------------------------------------------------+
| :ref:`udf.filter`             | Series or DataFrame    | Boolean                  | Only accepts UDFs in group by. Function is called for each group, and the group is removed from the result if the function returns ``False`` |
+-------------------------------+------------------------+--------------------------+----------------------------------------------------------------------------------------------------------------------------------------------+
| :ref:`udf.agg`                | Series or DataFrame    | Scalar or Series         | Aggregate and summarizes values, e.g., sum or custom reducer                                                                                 |
+-------------------------------+------------------------+--------------------------+----------------------------------------------------------------------------------------------------------------------------------------------+
| :ref:`udf.transform` (axis=0) | Column (Series)        | Column (Series)          | Same as :meth:`apply` with (axis=0), but it raises an exception if the function changes the shape of the data                                |
+-------------------------------+------------------------+--------------------------+----------------------------------------------------------------------------------------------------------------------------------------------+
| :ref:`udf.transform` (axis=1) | Row (Series)           | Row (Series)             | Same as :meth:`apply` with (axis=1), but it raises an exception if the function changes the shape of the data                                |
+-------------------------------+------------------------+--------------------------+----------------------------------------------------------------------------------------------------------------------------------------------+

When applying UDFs in pandas, it is essential to select the appropriate method based
on your specific task. Each method has its strengths and is designed for different use
cases. Understanding the purpose and behavior of each method will help you make informed
decisions, ensuring more efficient and maintainable code.

.. note::
    Some of these methods are can also be applied to groupby, resample, and various window objects.
    See :ref:`groupby`, :ref:`resample()<timeseries>`, :ref:`rolling()<window>`, :ref:`expanding()<window>`,
    and :ref:`ewm()<window>` for details.


.. _udf.map:

:meth:`Series.map` and :meth:`DataFrame.map`
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The :meth:`map` method is used specifically to apply element-wise UDFs. This means the function
will be called for each element in the ``Series`` or ``DataFrame``, with the individual value or
the cell as the function argument.

.. ipython:: python

    temperature_celsius = pd.DataFrame({
        "NYC": [14, 21, 23],
        "Los Angeles": [22, 28, 31],
    })

    def to_fahrenheit(value):
        return value * (9 / 5) + 32

    temperature_celsius.map(to_fahrenheit)

In this example, the function ``to_fahrenheit`` will be called 6 times, once for each value
in the ``DataFrame``. And the result of each call will be returned in the corresponding cell
of the resulting ``DataFrame``.

In general, ``map`` will be slow, as it will not make use of vectorization. Instead, a Python
function call for each value will be required, which will slow down things significantly if
working with medium or large data.

When to use: Use :meth:`map` for applying element-wise UDFs to DataFrames or Series.

.. _udf.apply:

:meth:`Series.apply` and :meth:`DataFrame.apply`
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The :meth:`apply` method allows you to apply UDFs for a whole column or row. This is different
from :meth:`map` in that the function will be called for each column (or row), not for each individual value.

.. ipython:: python

    temperature_celsius = pd.DataFrame({
        "NYC": [14, 21, 23],
        "Los Angeles": [22, 28, 31],
    })

    def to_fahrenheit(column):
        return column * (9 / 5) + 32

    temperature_celsius.apply(to_fahrenheit)

In the example, ``to_fahrenheit`` will be called only twice, as opposed to the 6 times with :meth:`map`.
This will be faster than using :meth:`map`, since the operations for each column are vectorized, and the
overhead of iterating over data in Python and calling Python functions is significantly reduced.

In some cases, the function may require all the data to be able to compute the result. So :meth:`apply`
is needed, since with :meth:`map` the function can only access one element at a time.

.. ipython:: python

    temperature = pd.DataFrame({
        "NYC": [14, 21, 23],
        "Los Angeles": [22, 28, 31],
    })

    def normalize(column):
        return column / column.mean()

    temperature.apply(normalize)

In the example, the ``normalize`` function needs to compute the mean of the whole column in order
to divide each element by it. So, we cannot call the function for each element, but we need the
function to receive the whole column.

:meth:`apply` can also execute function by row, by specifying ``axis=1``.

.. ipython:: python

    temperature = pd.DataFrame({
        "NYC": [14, 21, 23],
        "Los Angeles": [22, 28, 31],
    })

    def hotter(row):
        return row["Los Angeles"] - row["NYC"]

    temperature.apply(hotter, axis=1)

In the example, the function ``hotter`` will be called 3 times, once for each row. And each
call will receive the whole row as the argument, allowing computations that require more than
one value in the row.

``apply`` is also available for :meth:`SeriesGroupBy.apply`, :meth:`DataFrameGroupBy.apply`,
:meth:`Rolling.apply`, :meth:`Expanding.apply` and :meth:`Resampler.apply`. You can read more
about ``apply`` in groupby operations :ref:`groupby.apply`.

When to use: :meth:`apply` is suitable when no alternative vectorized method or UDF method is available,
but consider optimizing performance with vectorized operations wherever possible.

.. _udf.pipe:

:meth:`Series.pipe` and :meth:`DataFrame.pipe`
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The ``pipe`` method is similar to ``map`` and ``apply``, but the function receives the whole ``Series``
or ``DataFrame`` it is called on.

.. ipython:: python

    temperature = pd.DataFrame({
        "NYC": [14, 21, 23],
        "Los Angeles": [22, 28, 31],
    })

    def normalize(df):
        return df / df.mean().mean()

    temperature.pipe(normalize)

This is equivalent to calling the ``normalize`` function with the ``DataFrame`` as the parameter.

.. ipython:: python

    normalize(temperature)

The main advantage of using ``pipe`` is readability. It allows method chaining and clearer code when
calling multiple functions.

.. ipython:: python

    temperature_celsius = pd.DataFrame({
        "NYC": [14, 21, 23],
        "Los Angeles": [22, 28, 31],
    })

    def multiply_by_9(value):
        return value * 9

    def divide_by_5(value):
        return value / 5

    def add_32(value):
        return value + 32

    # Without `pipe`:
    fahrenheit = add_32(divide_by_5(multiply_by_9(temperature_celsius)))

    # With `pipe`:
    fahrenheit = (temperature_celsius.pipe(multiply_by_9)
                                     .pipe(divide_by_5)
                                     .pipe(add_32))

``pipe`` is also available for :meth:`SeriesGroupBy.pipe`, :meth:`DataFrameGroupBy.pipe` and
:meth:`Resampler.pipe`. You can read more about ``pipe`` in groupby operations in :ref:`groupby.pipe`.

When to use: Use :meth:`pipe` when you need to create a pipeline of operations and want to keep the code readable and maintainable.

.. _udf.filter:

:meth:`Series.filter` and :meth:`DataFrame.filter`
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The ``filter`` method is used to select a subset of rows that match certain criteria.
:meth:`Series.filter` and :meth:`DataFrame.filter` do not support user defined functions,
but :meth:`SeriesGroupBy.filter` and :meth:`DataFrameGroupBy.filter` do. You can read more
about ``filter`` in groupby operations in :ref:`groupby.filter`.

.. _udf.agg:

:meth:`Series.agg` and :meth:`DataFrame.agg`
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The ``agg`` method is used to aggregate a set of data points into a single one.
The most common aggregation functions such as ``min``, ``max``, ``mean``, ``sum``, etc.
are already implemented in pandas. ``agg`` allows to implement other custom aggregate
functions.

.. ipython:: python

    temperature = pd.DataFrame({
        "NYC": [14, 21, 23],
        "Los Angeles": [22, 28, 31],
    })

    def highest_jump(column):
        return column.pct_change().max()

    temperature.agg(highest_jump)


When to use: Use :meth:`agg` for performing custom aggregations, where the operation returns
a scalar value on each input.

.. _udf.transform:

:meth:`Series.transform` and :meth:`DataFrame.transform`
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The ``transform``` method is similar to an aggregation, with the difference that the result is broadcasted
to the original data.

.. ipython:: python

    temperature = pd.DataFrame({
        "NYC": [14, 21, 23],
        "Los Angeles": [22, 28, 31]},
        index=pd.date_range("2000-01-01", "2000-01-03"))

    def warm_up_all_days(column):
        return pd.Series(column.max(), index=column.index)

    temperature.transform(warm_up_all_days)

In the example, the ``warm_up_all_days`` function computes the ``max`` like an aggregation, but instead
of returning just the maximum value, it returns a ``DataFrame`` with the same shape as the original one
with the values of each day replaced by the maximum temperature of the city.

``transform`` is also available for :meth:`SeriesGroupBy.transform`, :meth:`DataFrameGroupBy.transform` and
:meth:`Resampler.transform`, where it's more common. You can read more about ``transform`` in groupby
operations in :ref:`groupby.transform`.

When to use: When you need to perform an aggregation that will be returned in the original structure of
the DataFrame.


Performance
-----------

While UDFs provide flexibility, their use is generally discouraged as they can introduce
performance issues, especially when written in pure Python. To improve efficiency,
consider using built-in ``NumPy`` or ``pandas`` functions instead of UDFs
for common operations.

.. note::
    If performance is critical, explore **vectorized operations** before resorting
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
via ``NumPy`` to process entire arrays at once. This approach avoids the overhead of looping
through rows in Python and making separate function calls for each row, which is slow and
inefficient. Additionally, ``NumPy`` arrays benefit from memory efficiency and CPU-level
optimizations, making vectorized operations the preferred choice whenever possible.


Improving Performance with UDFs
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

In scenarios where UDFs are necessary, there are still ways to mitigate their performance drawbacks.
One approach is to use **Numba**, a Just-In-Time (JIT) compiler that can significantly speed up numerical
Python code by compiling Python functions to optimized machine code at runtime.

By annotating your UDFs with ``@numba.jit``, you can achieve performance closer to vectorized operations,
especially for computationally heavy tasks.

.. note::
    You may also refer to the user guide on `Enhancing performance <https://pandas.pydata.org/pandas-docs/dev/user_guide/enhancingperf.html#numba-jit-compilation>`_
    for a more detailed guide to using **Numba**.

Using :meth:`DataFrame.pipe` for Composable Logic
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Another useful pattern for improving readability and composability, especially when mixing
vectorized logic with UDFs, is to use the :meth:`DataFrame.pipe` method.

:meth:`DataFrame.pipe` doesn't improve performance directly, but it enables cleaner
method chaining by passing the entire object into a function. This is especially helpful
when chaining custom transformations:

.. code-block:: python

    def add_ratio_column(df):
        df["ratio"] = 100 * (df["one"] / df["two"])
        return df

    df = (
        df
        .query("one > 0")
        .pipe(add_ratio_column)
        .dropna()
    )

This is functionally equivalent to calling ``add_ratio_column(df)``, but keeps your code
clean and composable. The function you pass to :meth:`DataFrame.pipe` can use vectorized operations,
row-wise UDFs, or any other logic; :meth:`DataFrame.pipe` is agnostic.

.. note::
    While :meth:`DataFrame.pipe` does not improve performance on its own,
    it promotes clean, modular design and allows both vectorized and UDF-based logic
    to be composed in method chains.
