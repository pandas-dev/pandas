"""
Implementation of pandas.Series.str and its interface.

* strings.accessor.StringMethods : Accessor for Series.str

Most methods on the StringMethods accessor follow the pattern:

    1. extract the array from the series (or index)
    2. Call that array's implementation of the string method
    3. Wrap the result (in a Series, index, or DataFrame)

To avoid namespace clashes and pollution,
these are prefixed with `_str_`. So ``Series.str.upper()`` calls
``Series.array._str_upper()``. The interface isn't currently public
to other string extension arrays.
"""
# Pandas current implementation is in ObjectStringArrayMixin. This is designed
# to work on object-dtype ndarrays.
#
#  - ObjectStringArrayMixin
#     - StringArray
#     - NumpyExtensionArray
#     - Categorical
#     - ArrowStringArray
