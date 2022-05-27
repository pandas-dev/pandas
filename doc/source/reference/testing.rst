{{ header }}

.. _api.testing:

=======
Testing
=======
.. currentmodule:: pandas

.. _api.general.testing:

Assertion functions
-------------------
.. autosummary::
   :toctree: api/

   testing.assert_frame_equal
   testing.assert_series_equal
   testing.assert_index_equal
   testing.assert_extension_array_equal

Exceptions and warnings
-----------------------
.. autosummary::
   :toctree: api/

   errors.AbstractMethodError
   errors.AccessorRegistrationWarning
   errors.DataError
   errors.DtypeWarning
   errors.DuplicateLabelError
   errors.EmptyDataError
   errors.InvalidIndexError
   errors.IntCastingNaNError
   errors.MergeError
   errors.NullFrequencyError
   errors.NumbaUtilError
   errors.OptionError
   errors.OutOfBoundsDatetime
   errors.OutOfBoundsTimedelta
   errors.ParserError
   errors.ParserWarning
   errors.PerformanceWarning
   errors.SettingWithCopyError
   errors.SpecificationError
   errors.UnsortedIndexError
   errors.UnsupportedFunctionCall

Bug report function
-------------------
.. autosummary::
   :toctree: api/

   show_versions

Test suite runner
-----------------
.. autosummary::
   :toctree: api/

   test
