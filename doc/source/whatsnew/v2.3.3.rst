.. _whatsnew_233:

What's new in 2.3.3 (September 29, 2025)
----------------------------------------

These are the changes in pandas 2.3.3. See :ref:`release` for a full changelog
including other versions of pandas.

{{ header }}

.. _whatsnew_233.py14_compat:

Pandas 2.3.3 is now compatible with Python 3.14
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Pandas 2.3.3 is the first version of pandas that is generally compatible with the upcoming
Python 3.14, and both wheels for free-threaded and normal Python 3.14 will be uploaded for
this release.

As usual please report any bugs discovered to our `issue tracker <https://github.com/pandas-dev/pandas/issues/new/choose>`_

.. ---------------------------------------------------------------------------
.. _whatsnew_233.string_fixes:

Improvements and fixes for the StringDtype
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Most changes in this release are related to :class:`StringDtype` which will
become the default string dtype in pandas 3.0. See
:ref:`whatsnew_230.upcoming_changes` for more details.

.. _whatsnew_233.string_fixes.improvements:

Improvements
^^^^^^^^^^^^
- Update :meth:`DataFrame.select_dtypes` to keep selecting ``str`` columns when
  specifying ``include=["object"]`` for backwards compatibility. In a future
  release, this will be deprecated and code for pandas 3+ should be updated to
  do ``include=["str"]`` (:issue:`61916`)
- Support the ``/`` operation between a ``pathlib.Path`` object and a :class:`StringDtype`
  Series, similarly as it works for object-dtype Series (:issue:`61940`)

.. _whatsnew_233.string_fixes.bugs:

Bug fixes
^^^^^^^^^
- Fix bug in :meth:`Series.str.replace` using named capture groups (e.g., ``\g<name>``) with the Arrow-backed dtype would raise an error (:issue:`57636`)
- Fix regression in :meth:`Series.str.contains`, :meth:`~Series.str.match` and :meth:`~Series.str.fullmatch`
  with a compiled regex and custom flags (:issue:`62240`)
- Fix :meth:`Series.str.match` and :meth:`~Series.str.fullmatch` not matching patterns with groups correctly for the Arrow-backed string dtype (:issue:`61072`)
- Fix bug in :meth:`~DataFrame.groupby` with ``sum()`` and unobserved categories resulting in ``0`` instead of the empty string ``""`` (:issue:`61909`)
- Fix :meth:`Series.str.isdigit` to correctly recognize unicode superscript
  characters as digits for :class:`StringDtype` backed by PyArrow (:issue:`61466`)
- Fix comparing a :class:`StringDtype` Series with mixed objects raising an error (:issue:`60228`)
- Fix error being raised when using a numpy ufunc with a Python-backed string array (:issue:`40800`)

Other changes
~~~~~~~~~~~~~

- The deprecation of using :meth:`Series.resample` and :meth:`DataFrame.resample`
  with a :class:`PeriodIndex` (and the 'convention' keyword) has been undone.
  Resampling with a :class:`PeriodIndex` is supported again, but a subset of
  methods that return incorrect results will raise an error in pandas 3.0 (:issue:`57033`)

Other bug fixes
~~~~~~~~~~~~~~~~

- Fix memory leak in :meth:`DataFrame.to_json` with datetime columns (:issue:`62204`)
- Fixed regression in :meth:`DataFrame.from_records` not initializing subclasses properly (:issue:`57008`)
- The :meth:`DataFrame.iloc` now works correctly with ``copy_on_write`` option when assigning values after subsetting the columns of a homogeneous DataFrame (:issue:`60309`)

.. ---------------------------------------------------------------------------
.. _whatsnew_233.contributors:

Contributors
~~~~~~~~~~~~

.. contributors:: v2.3.2..v2.3.3|HEAD
