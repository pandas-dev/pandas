.. _user_defined_functions:

{{ header }}

**************************************
Introduction to User Defined Functions
**************************************

In pandas, User Defined Functions (UDFs) provide a way to extend the library’s
functionality by allowing users to apply custom computations to their data. While
pandas comes with a set of built-in functions for data manipulation, UDFs offer
flexibility when built-in methods are not sufficient. These functions can be 
applied at different levels: element-wise, row-wise, column-wise, or group-wise,
depending on the method used.

Note: User Defined Functions will be abbreviated to UDFs throughout this guide.

Why Use UDFs?
-------------

Pandas is designed for high-performance data processing, but sometimes your specific
needs go beyond standard aggregation, transformation, or filtering. UDFs allow you to:
* Customize Computations: Implement logic tailored to your dataset, such as complex 
  transformations, domain-specific calculations, or conditional modifications.
* Improve Code Readability: Encapsulate logic into functions rather than writing long,
  complex expressions.
* Handle Complex Grouped Operations: Perform operations on grouped data that standard
  methods do not support.
* Extend pandas' Functionality: Apply external libraries or advanced calculations that 
  are not natively available.


Where Can UDFs Be Used?
-----------------------

UDFs can be applied across various pandas methods that work with both Series and DataFrames:

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