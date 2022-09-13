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
   errors.AttributeConflictWarning
   errors.CategoricalConversionWarning
   errors.ClosedFileError
   errors.CSSWarning
   errors.DatabaseError
   errors.DataError
   errors.DtypeWarning
   errors.DuplicateLabelError
   errors.EmptyDataError
   errors.IncompatibilityWarning
   errors.IndexingError
   errors.InvalidColumnName
   errors.InvalidComparison
   errors.InvalidIndexError
   errors.InvalidVersion
   errors.IntCastingNaNError
   errors.LossySetitemError
   errors.MergeError
   errors.NoBufferPresent
   errors.NullFrequencyError
   errors.NumbaUtilError
   errors.NumExprClobberingError
   errors.OptionError
   errors.OutOfBoundsDatetime
   errors.OutOfBoundsTimedelta
   errors.ParserError
   errors.ParserWarning
   errors.PerformanceWarning
   errors.PossibleDataLossError
   errors.PossiblePrecisionLoss
   errors.PyperclipException
   errors.PyperclipWindowsException
   errors.SettingWithCopyError
   errors.SettingWithCopyWarning
   errors.SpecificationError
   errors.UndefinedVariableError
   errors.UnsortedIndexError
   errors.UnsupportedFunctionCall
   errors.ValueLabelTypeMismatch

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
