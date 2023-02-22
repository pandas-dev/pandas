#!/bin/bash
#
# Run checks related to code quality.
#
# This script is intended for both the CI and to check locally that code standards are
# respected. We run doctests here (currently some files only), and we
# validate formatting error in docstrings.
#
# Usage:
#   $ ./ci/code_checks.sh               # run all checks
#   $ ./ci/code_checks.sh code          # checks on imported code
#   $ ./ci/code_checks.sh doctests      # run doctests
#   $ ./ci/code_checks.sh docstrings    # validate docstring errors
#   $ ./ci/code_checks.sh single-docs   # check single-page docs build warning-free
#   $ ./ci/code_checks.sh notebooks     # check execution of documentation notebooks

[[ -z "$1" || "$1" == "code" || "$1" == "doctests" || "$1" == "docstrings" || "$1" == "single-docs" || "$1" == "notebooks" ]] || \
    { echo "Unknown command $1. Usage: $0 [code|doctests|docstrings|single-docs|notebooks]"; exit 9999; }

BASE_DIR="$(dirname $0)/.."
RET=0
CHECK=$1

function invgrep {
    # grep with inverse exist status and formatting for azure-pipelines
    #
    # This function works exactly as grep, but with opposite exit status:
    # - 0 (success) when no patterns are found
    # - 1 (fail) when the patterns are found
    #
    # This is useful for the CI, as we want to fail if one of the patterns
    # that we want to avoid is found by grep.
    grep -n "$@" | sed "s/^/$INVGREP_PREPEND/" | sed "s/$/$INVGREP_APPEND/" ; EXIT_STATUS=${PIPESTATUS[0]}
    return $((! $EXIT_STATUS))
}

if [[ "$GITHUB_ACTIONS" == "true" ]]; then
    INVGREP_PREPEND="##[error]"
fi

### CODE ###
if [[ -z "$CHECK" || "$CHECK" == "code" ]]; then

    MSG='Check import. No warnings, and blocklist some optional dependencies' ; echo $MSG
    python -W error -c "
import sys
import pandas

blocklist = {'bs4', 'gcsfs', 'html5lib', 'http', 'ipython', 'jinja2', 'hypothesis',
             'lxml', 'matplotlib', 'openpyxl', 'py', 'pytest', 's3fs', 'scipy',
             'tables', 'urllib.request', 'xlrd', 'xlsxwriter'}

# GH#28227 for some of these check for top-level modules, while others are
#  more specific (e.g. urllib.request)
import_mods = set(m.split('.')[0] for m in sys.modules) | set(sys.modules)
mods = blocklist & import_mods
if mods:
    sys.stderr.write('err: pandas should not import: {}\n'.format(', '.join(mods)))
    sys.exit(len(mods))
    "
    RET=$(($RET + $?)) ; echo $MSG "DONE"

fi

### DOCTESTS ###
if [[ -z "$CHECK" || "$CHECK" == "doctests" ]]; then

    MSG='Doctests' ; echo $MSG
    # Ignore test_*.py files or else the unit tests will run
    python -m pytest --doctest-modules --ignore-glob="**/test_*.py" pandas
    RET=$(($RET + $?)) ; echo $MSG "DONE"

    MSG='Cython Doctests' ; echo $MSG
    python -m pytest --doctest-cython pandas/_libs
    RET=$(($RET + $?)) ; echo $MSG "DONE"

fi

### DOCSTRINGS ###
if [[ -z "$CHECK" || "$CHECK" == "docstrings" ]]; then

    MSG='Validate docstrings (EX04, GL01, GL02, GL03, GL04, GL05, GL06, GL07, GL09, GL10, PR03, PR04, PR05, PR06, PR08, PR09, PR10, RT01, RT02, RT04, RT05, SA02, SA03, SA04, SS01, SS02, SS03, SS04, SS05, SS06)' ; echo $MSG
    $BASE_DIR/scripts/validate_docstrings.py --format=actions --errors=EX04,GL01,GL02,GL03,GL04,GL05,GL06,GL07,GL09,GL10,PR03,PR04,PR05,PR06,PR08,PR09,PR10,RT01,RT02,RT04,RT05,SA02,SA03,SA04,SS01,SS02,SS03,SS04,SS05,SS06
    RET=$(($RET + $?)) ; echo $MSG "DONE"

    MSG='Partially validate docstrings (EX01)' ;  echo $MSG
    $BASE_DIR/scripts/validate_docstrings.py --format=actions --errors=EX01 --ignore_functions \
        pandas.Series.index \
        pandas.Series.hasnans \
        pandas.Series.to_list \
        pandas.Series.__iter__ \
        pandas.Series.keys \
        pandas.Series.item \
        pandas.Series.pipe \
        pandas.Series.kurt \
        pandas.Series.mean \
        pandas.Series.median \
        pandas.Series.mode \
        pandas.Series.sem \
        pandas.Series.skew \
        pandas.Series.kurtosis \
        pandas.Series.is_unique \
        pandas.Series.is_monotonic_increasing \
        pandas.Series.is_monotonic_decreasing \
        pandas.Series.backfill \
        pandas.Series.pad \
        pandas.Series.argsort \
        pandas.Series.reorder_levels \
        pandas.Series.ravel \
        pandas.Series.first_valid_index \
        pandas.Series.last_valid_index \
        pandas.Series.dt.date \
        pandas.Series.dt.time \
        pandas.Series.dt.timetz \
        pandas.Series.dt.dayofyear \
        pandas.Series.dt.day_of_year \
        pandas.Series.dt.quarter \
        pandas.Series.dt.daysinmonth \
        pandas.Series.dt.days_in_month \
        pandas.Series.dt.tz \
        pandas.Series.dt.end_time \
        pandas.Series.dt.days \
        pandas.Series.dt.seconds \
        pandas.Series.dt.microseconds \
        pandas.Series.dt.nanoseconds \
        pandas.Series.str.center \
        pandas.Series.str.decode \
        pandas.Series.str.encode \
        pandas.Series.str.find \
        pandas.Series.str.fullmatch \
        pandas.Series.str.index \
        pandas.Series.str.ljust \
        pandas.Series.str.match \
        pandas.Series.str.normalize \
        pandas.Series.str.rfind \
        pandas.Series.str.rindex \
        pandas.Series.str.rjust \
        pandas.Series.str.translate \
        pandas.Series.sparse \
        pandas.DataFrame.sparse \
        pandas.Series.cat.categories \
        pandas.Series.cat.ordered \
        pandas.Series.cat.codes \
        pandas.Series.cat.reorder_categories \
        pandas.Series.cat.set_categories \
        pandas.Series.cat.as_ordered \
        pandas.Series.cat.as_unordered \
        pandas.Series.sparse.fill_value \
        pandas.Flags \
        pandas.Series.attrs \
        pandas.Series.plot \
        pandas.Series.hist \
        pandas.Series.to_string \
        pandas.errors.AbstractMethodError \
        pandas.errors.AccessorRegistrationWarning \
        pandas.errors.AttributeConflictWarning \
        pandas.errors.DataError \
        pandas.errors.EmptyDataError \
        pandas.errors.IncompatibilityWarning \
        pandas.errors.InvalidComparison \
        pandas.errors.InvalidIndexError \
        pandas.errors.InvalidVersion \
        pandas.errors.IntCastingNaNError \
        pandas.errors.LossySetitemError \
        pandas.errors.MergeError \
        pandas.errors.NoBufferPresent \
        pandas.errors.NullFrequencyError \
        pandas.errors.NumbaUtilError \
        pandas.errors.OptionError \
        pandas.errors.OutOfBoundsDatetime \
        pandas.errors.OutOfBoundsTimedelta \
        pandas.errors.ParserError \
        pandas.errors.PerformanceWarning \
        pandas.errors.PyperclipException \
        pandas.errors.PyperclipWindowsException \
        pandas.errors.UnsortedIndexError \
        pandas.errors.UnsupportedFunctionCall \
        pandas.show_versions \
        pandas.test \
        pandas.NaT \
        pandas.Timestamp.as_unit \
        pandas.Timestamp.ctime \
        pandas.Timestamp.date \
        pandas.Timestamp.dst \
        pandas.Timestamp.isocalendar \
        pandas.Timestamp.isoweekday \
        pandas.Timestamp.strptime \
        pandas.Timestamp.time \
        pandas.Timestamp.timetuple \
        pandas.Timestamp.timetz \
        pandas.Timestamp.to_datetime64 \
        pandas.Timestamp.toordinal \
        pandas.Timestamp.tzname \
        pandas.Timestamp.utcoffset \
        pandas.Timestamp.utctimetuple \
        pandas.Timestamp.weekday \
        pandas.arrays.DatetimeArray \
        pandas.Timedelta.view \
        pandas.Timedelta.as_unit \
        pandas.Timedelta.ceil \
        pandas.Timedelta.floor \
        pandas.Timedelta.round \
        pandas.Timedelta.to_pytimedelta \
        pandas.Timedelta.to_timedelta64 \
        pandas.Timedelta.to_numpy \
        pandas.Timedelta.total_seconds \
        pandas.arrays.TimedeltaArray \
        pandas.Period.end_time \
        pandas.Period.freqstr \
        pandas.Period.is_leap_year \
        pandas.Period.month \
        pandas.Period.quarter \
        pandas.Period.year \
        pandas.Period.asfreq \
        pandas.Period.now \
        pandas.arrays.PeriodArray \
        pandas.Interval.closed \
        pandas.Interval.left \
        pandas.Interval.length \
        pandas.Interval.right \
        pandas.arrays.IntervalArray.left \
        pandas.arrays.IntervalArray.right \
        pandas.arrays.IntervalArray.closed \
        pandas.arrays.IntervalArray.mid \
        pandas.arrays.IntervalArray.length \
        pandas.arrays.IntervalArray.is_non_overlapping_monotonic \
        pandas.arrays.IntervalArray.from_arrays \
        pandas.arrays.IntervalArray.to_tuples \
        pandas.Int8Dtype \
        pandas.Int16Dtype \
        pandas.Int32Dtype \
        pandas.Int64Dtype \
        pandas.UInt8Dtype \
        pandas.UInt16Dtype \
        pandas.UInt32Dtype \
        pandas.UInt64Dtype \
        pandas.NA \
        pandas.Float32Dtype \
        pandas.Float64Dtype \
        pandas.CategoricalDtype.categories \
        pandas.CategoricalDtype.ordered \
        pandas.Categorical.dtype \
        pandas.Categorical.categories \
        pandas.Categorical.ordered \
        pandas.Categorical.codes \
        pandas.Categorical.__array__ \
        pandas.SparseDtype \
        pandas.DatetimeTZDtype.unit \
        pandas.DatetimeTZDtype.tz \
        pandas.PeriodDtype.freq \
        pandas.IntervalDtype.subtype \
        pandas_dtype \
        pandas.api.types.is_bool \
        pandas.api.types.is_complex \
        pandas.api.types.is_float \
        pandas.api.types.is_integer \
        pandas.api.types.pandas_dtype \
        pandas.read_clipboard \
        pandas.ExcelFile \
        pandas.ExcelFile.parse \
        pandas.DataFrame.to_html \
        pandas.io.formats.style.Styler.to_html \
        pandas.HDFStore.put \
        pandas.HDFStore.append \
        pandas.HDFStore.get \
        pandas.HDFStore.select \
        pandas.HDFStore.info \
        pandas.HDFStore.keys \
        pandas.HDFStore.groups \
        pandas.HDFStore.walk \
        pandas.read_feather \
        pandas.DataFrame.to_feather \
        pandas.read_parquet \
        pandas.read_orc \
        pandas.read_sas \
        pandas.read_spss \
        pandas.read_sql_query \
        pandas.read_gbq \
        pandas.io.stata.StataReader.data_label \
        pandas.io.stata.StataReader.value_labels \
        pandas.io.stata.StataReader.variable_labels \
        pandas.io.stata.StataWriter.write_file \
        pandas.core.resample.Resampler.__iter__ \
        pandas.core.resample.Resampler.groups \
        pandas.core.resample.Resampler.indices \
        pandas.core.resample.Resampler.get_group \
        pandas.core.resample.Resampler.ffill \
        pandas.core.resample.Resampler.asfreq \
        pandas.core.resample.Resampler.count \
        pandas.core.resample.Resampler.nunique \
        pandas.core.resample.Resampler.max \
        pandas.core.resample.Resampler.mean \
        pandas.core.resample.Resampler.median \
        pandas.core.resample.Resampler.min \
        pandas.core.resample.Resampler.ohlc \
        pandas.core.resample.Resampler.prod \
        pandas.core.resample.Resampler.size \
        pandas.core.resample.Resampler.sem \
        pandas.core.resample.Resampler.std \
        pandas.core.resample.Resampler.sum \
        pandas.core.resample.Resampler.var \
        pandas.core.resample.Resampler.quantile \
        pandas.describe_option \
        pandas.reset_option \
        pandas.get_option \
        pandas.set_option \
        pandas.plotting.deregister_matplotlib_converters \
        pandas.plotting.plot_params \
        pandas.plotting.register_matplotlib_converters \
        pandas.plotting.table \
        pandas.util.hash_array \
        pandas.util.hash_pandas_object \
        pandas_object \
        pandas.api.interchange.from_dataframe \
        pandas.Index.values \
        pandas.Index.hasnans \
        pandas.Index.dtype \
        pandas.Index.inferred_type \
        pandas.Index.shape \
        pandas.Index.name \
        pandas.Index.nbytes \
        pandas.Index.ndim \
        pandas.Index.size \
        pandas.Index.T \
        pandas.Index.memory_usage \
        pandas.Index.copy \
        pandas.Index.drop \
        pandas.Index.identical \
        pandas.Index.insert \
        pandas.Index.is_ \
        pandas.Index.take \
        pandas.Index.putmask \
        pandas.Index.unique \
        pandas.Index.fillna \
        pandas.Index.dropna \
        pandas.Index.astype \
        pandas.Index.item \
        pandas.Index.map \
        pandas.Index.ravel \
        pandas.Index.to_list \
        pandas.Index.append \
        pandas.Index.join \
        pandas.Index.asof_locs \
        pandas.Index.get_slice_bound \
        pandas.RangeIndex \
        pandas.RangeIndex.start \
        pandas.RangeIndex.stop \
        pandas.RangeIndex.step \
        pandas.RangeIndex.from_range \
        pandas.CategoricalIndex.codes \
        pandas.CategoricalIndex.categories \
        pandas.CategoricalIndex.ordered \
        pandas.CategoricalIndex.reorder_categories \
        pandas.CategoricalIndex.set_categories \
        pandas.CategoricalIndex.as_ordered \
        pandas.CategoricalIndex.as_unordered \
        pandas.CategoricalIndex.equals \
        pandas.IntervalIndex.closed \
        pandas.IntervalIndex.values \
        pandas.IntervalIndex.is_non_overlapping_monotonic \
        pandas.IntervalIndex.to_tuples \
        pandas.MultiIndex.dtypes \
        pandas.MultiIndex.drop \
        pandas.DatetimeIndex \
        pandas.DatetimeIndex.date \
        pandas.DatetimeIndex.time \
        pandas.DatetimeIndex.timetz \
        pandas.DatetimeIndex.dayofyear \
        pandas.DatetimeIndex.day_of_year \
        pandas.DatetimeIndex.quarter \
        pandas.DatetimeIndex.tz \
        pandas.DatetimeIndex.freqstr \
        pandas.DatetimeIndex.inferred_freq \
        pandas.DatetimeIndex.indexer_at_time \
        pandas.DatetimeIndex.indexer_between_time \
        pandas.DatetimeIndex.snap \
        pandas.DatetimeIndex.as_unit \
        pandas.DatetimeIndex.to_pydatetime \
        pandas.DatetimeIndex.to_series \
        pandas.DatetimeIndex.mean \
        pandas.DatetimeIndex.std \
        pandas.TimedeltaIndex \
        pandas.TimedeltaIndex.days \
        pandas.TimedeltaIndex.seconds \
        pandas.TimedeltaIndex.microseconds \
        pandas.TimedeltaIndex.nanoseconds \
        pandas.TimedeltaIndex.components \
        pandas.TimedeltaIndex.inferred_freq \
        pandas.TimedeltaIndex.as_unit \
        pandas.TimedeltaIndex.to_pytimedelta \
        pandas.TimedeltaIndex.mean \
        pandas.PeriodIndex.day \
        pandas.PeriodIndex.dayofweek \
        pandas.PeriodIndex.day_of_week \
        pandas.PeriodIndex.dayofyear \
        pandas.PeriodIndex.day_of_year \
        pandas.PeriodIndex.days_in_month \
        pandas.PeriodIndex.daysinmonth \
        pandas.PeriodIndex.end_time \
        pandas.PeriodIndex.freqstr \
        pandas.PeriodIndex.hour \
        pandas.PeriodIndex.is_leap_year \
        pandas.PeriodIndex.minute \
        pandas.PeriodIndex.month \
        pandas.PeriodIndex.quarter \
        pandas.PeriodIndex.second \
        pandas.PeriodIndex.week \
        pandas.PeriodIndex.weekday \
        pandas.PeriodIndex.weekofyear \
        pandas.PeriodIndex.year \
        pandas.PeriodIndex.to_timestamp \
        pandas.core.window.rolling.Rolling.max \
        pandas.core.window.rolling.Rolling.cov \
        pandas.core.window.rolling.Rolling.skew \
        pandas.core.window.rolling.Rolling.apply \
        pandas.core.window.rolling.Window.mean \
        pandas.core.window.rolling.Window.sum \
        pandas.core.window.rolling.Window.var \
        pandas.core.window.rolling.Window.std \
        pandas.core.window.expanding.Expanding.count \
        pandas.core.window.expanding.Expanding.sum \
        pandas.core.window.expanding.Expanding.mean \
        pandas.core.window.expanding.Expanding.median \
        pandas.core.window.expanding.Expanding.min \
        pandas.core.window.expanding.Expanding.max \
        pandas.core.window.expanding.Expanding.corr \
        pandas.core.window.expanding.Expanding.cov \
        pandas.core.window.expanding.Expanding.skew \
        pandas.core.window.expanding.Expanding.apply \
        pandas.core.window.expanding.Expanding.quantile \
        pandas.core.window.ewm.ExponentialMovingWindow.mean \
        pandas.core.window.ewm.ExponentialMovingWindow.sum \
        pandas.core.window.ewm.ExponentialMovingWindow.std \
        pandas.core.window.ewm.ExponentialMovingWindow.var \
        pandas.core.window.ewm.ExponentialMovingWindow.corr \
        pandas.core.window.ewm.ExponentialMovingWindow.cov \
        pandas.api.indexers.BaseIndexer \
        pandas.api.indexers.VariableOffsetWindowIndexer \
        pandas.core.groupby.DataFrameGroupBy.__iter__ \
        pandas.core.groupby.SeriesGroupBy.__iter__ \
        pandas.core.groupby.DataFrameGroupBy.groups \
        pandas.core.groupby.SeriesGroupBy.groups \
        pandas.core.groupby.DataFrameGroupBy.indices \
        pandas.core.groupby.SeriesGroupBy.indices \
        pandas.core.groupby.DataFrameGroupBy.get_group \
        pandas.core.groupby.SeriesGroupBy.get_group \
        pandas.core.groupby.DataFrameGroupBy.all \
        pandas.core.groupby.DataFrameGroupBy.any \
        pandas.core.groupby.DataFrameGroupBy.bfill \
        pandas.core.groupby.DataFrameGroupBy.count \
        pandas.core.groupby.DataFrameGroupBy.cummax \
        pandas.core.groupby.DataFrameGroupBy.cummin \
        pandas.core.groupby.DataFrameGroupBy.cumprod \
        pandas.core.groupby.DataFrameGroupBy.cumsum \
        pandas.core.groupby.DataFrameGroupBy.diff \
        pandas.core.groupby.DataFrameGroupBy.ffill \
        pandas.core.groupby.DataFrameGroupBy.max \
        pandas.core.groupby.DataFrameGroupBy.median \
        pandas.core.groupby.DataFrameGroupBy.min \
        pandas.core.groupby.DataFrameGroupBy.ohlc \
        pandas.core.groupby.DataFrameGroupBy.pct_change \
        pandas.core.groupby.DataFrameGroupBy.prod \
        pandas.core.groupby.DataFrameGroupBy.sem \
        pandas.core.groupby.DataFrameGroupBy.shift \
        pandas.core.groupby.DataFrameGroupBy.size \
        pandas.core.groupby.DataFrameGroupBy.skew \
        pandas.core.groupby.DataFrameGroupBy.std \
        pandas.core.groupby.DataFrameGroupBy.sum \
        pandas.core.groupby.DataFrameGroupBy.var \
        pandas.core.groupby.SeriesGroupBy.all \
        pandas.core.groupby.SeriesGroupBy.any \
        pandas.core.groupby.SeriesGroupBy.bfill \
        pandas.core.groupby.SeriesGroupBy.count \
        pandas.core.groupby.SeriesGroupBy.cummax \
        pandas.core.groupby.SeriesGroupBy.cummin \
        pandas.core.groupby.SeriesGroupBy.cumprod \
        pandas.core.groupby.SeriesGroupBy.cumsum \
        pandas.core.groupby.SeriesGroupBy.diff \
        pandas.core.groupby.SeriesGroupBy.ffill \
        pandas.core.groupby.SeriesGroupBy.is_monotonic_increasing \
        pandas.core.groupby.SeriesGroupBy.is_monotonic_decreasing \
        pandas.core.groupby.SeriesGroupBy.max \
        pandas.core.groupby.SeriesGroupBy.median \
        pandas.core.groupby.SeriesGroupBy.min \
        pandas.core.groupby.SeriesGroupBy.nunique \
        pandas.core.groupby.SeriesGroupBy.ohlc \
        pandas.core.groupby.SeriesGroupBy.pct_change \
        pandas.core.groupby.SeriesGroupBy.prod \
        pandas.core.groupby.SeriesGroupBy.sem \
        pandas.core.groupby.SeriesGroupBy.shift \
        pandas.core.groupby.SeriesGroupBy.size \
        pandas.core.groupby.SeriesGroupBy.skew \
        pandas.core.groupby.SeriesGroupBy.std \
        pandas.core.groupby.SeriesGroupBy.sum \
        pandas.core.groupby.SeriesGroupBy.var \
        pandas.core.groupby.SeriesGroupBy.hist \
        pandas.core.groupby.DataFrameGroupBy.plot \
        pandas.core.groupby.SeriesGroupBy.plot \
        pandas.io.formats.style.Styler \
        pandas.io.formats.style.Styler.from_custom_template \
        pandas.io.formats.style.Styler.set_caption \
        pandas.io.formats.style.Styler.set_sticky \
        pandas.io.formats.style.Styler.set_uuid \
        pandas.io.formats.style.Styler.clear \
        pandas.io.formats.style.Styler.highlight_null \
        pandas.io.formats.style.Styler.highlight_max \
        pandas.io.formats.style.Styler.highlight_min \
        pandas.io.formats.style.Styler.bar \
        pandas.io.formats.style.Styler.to_string \
        pandas.api.extensions.ExtensionDtype \
        pandas.api.extensions.ExtensionArray \
        pandas.arrays.PandasArray \
        pandas.api.extensions.ExtensionArray._accumulate \
        pandas.api.extensions.ExtensionArray._concat_same_type \
        pandas.api.extensions.ExtensionArray._formatter \
        pandas.api.extensions.ExtensionArray._from_factorized \
        pandas.api.extensions.ExtensionArray._from_sequence \
        pandas.api.extensions.ExtensionArray._from_sequence_of_strings \
        pandas.api.extensions.ExtensionArray._reduce \
        pandas.api.extensions.ExtensionArray._values_for_argsort \
        pandas.api.extensions.ExtensionArray._values_for_factorize \
        pandas.api.extensions.ExtensionArray.argsort \
        pandas.api.extensions.ExtensionArray.astype \
        pandas.api.extensions.ExtensionArray.copy \
        pandas.api.extensions.ExtensionArray.view \
        pandas.api.extensions.ExtensionArray.dropna \
        pandas.api.extensions.ExtensionArray.equals \
        pandas.api.extensions.ExtensionArray.factorize \
        pandas.api.extensions.ExtensionArray.fillna \
        pandas.api.extensions.ExtensionArray.insert \
        pandas.api.extensions.ExtensionArray.isin \
        pandas.api.extensions.ExtensionArray.isna \
        pandas.api.extensions.ExtensionArray.ravel \
        pandas.api.extensions.ExtensionArray.searchsorted \
        pandas.api.extensions.ExtensionArray.shift \
        pandas.api.extensions.ExtensionArray.unique \
        pandas.api.extensions.ExtensionArray.dtype \
        pandas.api.extensions.ExtensionArray.nbytes \
        pandas.api.extensions.ExtensionArray.ndim \
        pandas.api.extensions.ExtensionArray.shape \
        pandas.api.extensions.ExtensionArray.tolist \
        pandas.DataFrame.index \
        pandas.DataFrame.columns \
        pandas.DataFrame.__iter__ \
        pandas.DataFrame.keys \
        pandas.DataFrame.iterrows \
        pandas.DataFrame.pipe \
        pandas.DataFrame.kurt \
        pandas.DataFrame.kurtosis \
        pandas.DataFrame.mean \
        pandas.DataFrame.median \
        pandas.DataFrame.sem \
        pandas.DataFrame.skew \
        pandas.DataFrame.backfill \
        pandas.DataFrame.pad \
        pandas.DataFrame.swapaxes \
        pandas.DataFrame.first_valid_index \
        pandas.DataFrame.last_valid_index \
        pandas.DataFrame.attrs \
        pandas.DataFrame.plot \
        pandas.DataFrame.sparse.density \
        pandas.DataFrame.sparse.to_coo \
        pandas.DataFrame.to_gbq \
        pandas.DataFrame.style \
        pandas.DataFrame.__dataframe__
    RET=$(($RET + $?)) ; echo $MSG "DONE"

    MSG='Partially validate docstrings (EX02)' ;  echo $MSG
    $BASE_DIR/scripts/validate_docstrings.py --format=actions --errors=EX02 --ignore_functions \
        pandas.DataFrame.plot.line \
        pandas.Series.plot.line \
        pandas.api.types.infer_dtype \
        pandas.api.types.is_datetime64_any_dtype \
        pandas.api.types.is_datetime64_ns_dtype \
        pandas.api.types.is_datetime64tz_dtype \
        pandas.api.types.is_integer_dtype \
        pandas.api.types.is_interval_dtype \
        pandas.api.types.is_period_dtype \
        pandas.api.types.is_signed_integer_dtype \
        pandas.api.types.is_sparse \
        pandas.api.types.is_string_dtype \
        pandas.api.types.is_unsigned_integer_dtype \
        pandas.plotting.andrews_curves \
        pandas.plotting.autocorrelation_plot \
        pandas.plotting.lag_plot \
        pandas.plotting.parallel_coordinates \
        pandas.plotting.radviz \
        pandas.tseries.frequencies.to_offset
    RET=$(($RET + $?)) ; echo $MSG "DONE"

fi

### DOCUMENTATION NOTEBOOKS ###
if [[ -z "$CHECK" || "$CHECK" == "notebooks" ]]; then

    MSG='Notebooks' ; echo $MSG
    jupyter nbconvert --execute $(find doc/source -name '*.ipynb') --to notebook
    RET=$(($RET + $?)) ; echo $MSG "DONE"

fi

### SINGLE-PAGE DOCS ###
if [[ -z "$CHECK" || "$CHECK" == "single-docs" ]]; then
    python doc/make.py --warnings-are-errors --single pandas.Series.value_counts
    python doc/make.py --warnings-are-errors --single pandas.Series.str.split
    python doc/make.py clean
fi

exit $RET
