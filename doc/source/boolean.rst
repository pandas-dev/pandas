.. currentmodule:: pandas

.. _boolean:

**************************
Nullable Boolean Data Type
**************************

.. versionadded:: 1.0.0

.. _boolean.klean:

Kleene Logic
------------

:class:`arrays.BooleanArray` implements Kleene logic (sometime called three-value logic) for
logical operations like ``&`` (and), ``|`` (or) and ``^`` (exclusive-or).

Here's a table for ``and``.

==========  ===========  ============
left value  right value  output value
==========  ===========  ============
True        True         True
True        False        False
True        NA           NA
False       False        False
False       NA           False
NA          NA           NA
==========  ===========  ============


And for ``or``

==========  ===========  ============
left value  right value  output value
==========  ===========  ============
True        True         True
True        False        True
True        NA           True
False       False        False
False       NA           NA
NA          NA           NA
==========  ===========  ============

And for ``xor``

==========  ===========  ============
left value  right value  output value
==========  ===========  ============
True        True         False
True        False        True
True        NA           NA
False       False        False
False       NA           NA
NA          NA           NA
==========  ===========  ============

When an ``NA`` is present in an operation, the output value is ``NA`` only if
the result cannot be determined soley based on the other input. For example,
``True | NA`` is ``True``, because both ``True | True`` and ``True | False``
are ``True``. In that case, we don't actually need to consider the value
of the ``NA``.

On the other hand, ``True & NA`` is ``NA``. The result depends on whether
the ``NA`` really is ``True`` or ``False``, since ``True & True`` is ``True``,
but ``True & False`` is ``False``, so we can't determine the output.


This differs from how ``np.nan`` behaves in logical operations. Pandas treated
``np.nan`` is *always false in the output*.

In ``or``

.. ipython:: python

   pd.Series([True, False, np.nan], dtype="object") | True
   pd.Series([True, False, np.nan], dtype="boolean") | True

In ``and``

   pd.Series([True, False, np.nan], dtype="object") & True
   pd.Series([True, False, np.nan], dtype="boolean") & True
