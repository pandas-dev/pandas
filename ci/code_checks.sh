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

set -uo pipefail

if [[ -v 1 ]]; then
    CHECK=$1
else
    # script will fail if it uses an unset variable (i.e. $1 is not provided)
    CHECK=""
fi

[[ -z "$CHECK" || "$CHECK" == "code" || "$CHECK" == "doctests" || "$CHECK" == "docstrings" || "$CHECK" == "single-docs" || "$CHECK" == "notebooks" ]] || \
    { echo "Unknown command $1. Usage: $0 [code|doctests|docstrings|single-docs|notebooks]"; exit 1; }

BASE_DIR="$(dirname "$0")/.."
RET=0

### CODE ###
if [[ -z "$CHECK" || "$CHECK" == "code" ]]; then

    MSG='Check import. No warnings, and blocklist some optional dependencies' ; echo "$MSG"
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
    RET=$(($RET + $?)) ; echo "$MSG" "DONE"

fi

### DOCTESTS ###
if [[ -z "$CHECK" || "$CHECK" == "doctests" ]]; then

    MSG='Python and Cython Doctests' ; echo "$MSG"
    # Using future.python_scalars=True avoids NumPy scalar reprs in docstrings.
    PANDAS_FUTURE_PYTHON_SCALARS="1" python -c 'import pandas as pd; pd.test(run_doctests=True)'
    RET=$(($RET + $?)) ; echo "$MSG" "DONE"

fi

### DOCSTRINGS ###
if [[ -z "$CHECK" || "$CHECK" == "docstrings" ]]; then

    MSG='Validate Docstrings' ; echo "$MSG"
    "$BASE_DIR"/scripts/validate_docstrings.py \
        --format=actions \
        -i "pandas.tseries.frequencies.to_offset ES01" \
        -i "pandas.Series.transform ES01" \
        -i "pandas.Series.rolling ES01" \
        -i "pandas.Series.count ES01" \
        -i "pandas.Series.mean ES01" \
        -i "pandas.Series.median ES01" \
        -i "pandas.Series.nlargest ES01" \
        -i "pandas.Series.nsmallest ES01" \
        -i "pandas.Series.quantile ES01" \
        -i "pandas.Series.is_unique ES01" \
        -i "pandas.Series.is_monotonic_increasing ES01" \
        -i "pandas.Series.is_monotonic_decreasing ES01" \
        -i "pandas.Series.drop_duplicates ES01" \
        -i "pandas.Series.rename_axis ES01" \
        -i "pandas.Series.unstack ES01" \
        -i "pandas.Series.explode ES01" \
        -i "pandas.Series.compare ES01" \
        -i "pandas.Series.at_time ES01" \
        -i "pandas.Series.cat ES01" \
        -i "pandas.Series.sparse ES01" \
        -i "pandas.Series.cat.codes ES01" \
        -i "pandas.Series.sparse.density ES01" \
        -i "pandas.Series.list.flatten ES01" \
        -i "pandas.Series.list.len ES01" \
        -i "pandas.Series.list.__getitem__ ES01" \
        -i "pandas.Series.struct.dtypes ES01" \
        -i "pandas.Series.struct.field ES01" \
        -i "pandas.Series.struct.explode ES01" \
        -i "pandas.Series.hist ES01" \
        -i "pandas.Series.to_pickle ES01" \
        -i "pandas.Series.to_csv ES01" \
        -i "pandas.Series.to_dict ES01" \
        -i "pandas.Series.to_frame ES01" \
        -i "pandas.Series.to_string ES01" \
        -i "pandas.Series.to_markdown ES01" \
        -i "pandas.testing.assert_series_equal ES01" \
        -i "pandas.testing.assert_index_equal ES01" \
        -i "pandas.Timestamp.is_month_end ES01" \
        -i "pandas.Timestamp.is_month_start ES01" \
        -i "pandas.Timestamp.is_quarter_end ES01" \
        -i "pandas.Timestamp.is_quarter_start ES01" \
        -i "pandas.Timestamp.is_year_end ES01" \
        -i "pandas.Timestamp.is_year_start ES01" \
        -i "pandas.Timestamp.microsecond ES01" \
        -i "pandas.Timestamp.minute ES01" \
        -i "pandas.Timestamp.month ES01" \
        -i "pandas.Timestamp.nanosecond ES01" \
        -i "pandas.Timestamp.second ES01" \
        -i "pandas.Timestamp.value ES01" \
        -i "pandas.Timestamp.week ES01" \
        -i "pandas.Timestamp.weekofyear ES01" \
        -i "pandas.Timestamp.year ES01" \
        -i "pandas.Timestamp.as_unit ES01" \
        -i "pandas.Timestamp.ceil ES01" \
        -i "pandas.Timestamp.day_name ES01" \
        -i "pandas.Timestamp.floor ES01" \
        -i "pandas.Timestamp.isocalendar ES01" \
        -i "pandas.Timestamp.strftime ES01" \
        -i "pandas.Timestamp.utcnow ES01" \
        -i "pandas.Timedelta.nanoseconds ES01" \
        -i "pandas.Timedelta.as_unit ES01" \
        -i "pandas.Timedelta.ceil ES01" \
        -i "pandas.Timedelta.floor ES01" \
        -i "pandas.Timedelta.round ES01" \
        -i "pandas.Interval ES01" \
        -i "pandas.Interval.left ES01" \
        -i "pandas.Interval.length ES01" \
        -i "pandas.Interval.mid ES01" \
        -i "pandas.Interval.right ES01" \
        -i "pandas.CategoricalDtype.categories ES01" \
        -i "pandas.CategoricalDtype.ordered ES01" \
        -i "pandas.DatetimeTZDtype.unit ES01" \
        -i "pandas.DatetimeTZDtype.tz ES01" \
        -i "pandas.IntervalDtype.subtype ES01" \
        -i "pandas.api.types.pandas_dtype ES01" \
        -i "pandas.api.types.is_any_real_numeric_dtype ES01" \
        -i "pandas.api.types.is_complex_dtype ES01" \
        -i "pandas.api.types.is_datetime64_any_dtype ES01" \
        -i "pandas.api.types.is_datetime64_dtype ES01" \
        -i "pandas.api.types.is_datetime64_ns_dtype ES01" \
        -i "pandas.api.types.is_dtype_equal ES01" \
        -i "pandas.api.types.is_numeric_dtype ES01" \
        -i "pandas.api.types.is_timedelta64_dtype ES01" \
        -i "pandas.api.types.is_dict_like ES01" \
        -i "pandas.api.types.is_named_tuple ES01" \
        -i "pandas.api.types.is_bool ES01" \
        -i "pandas.api.types.is_complex ES01" \
        -i "pandas.api.types.is_re ES01" \
        -i "pandas.api.types.is_re_compilable ES01" \
        -i "pandas.api.types.is_scalar ES01" \
        -i "pandas.DataFrame.to_pickle ES01" \
        -i "pandas.DataFrame.to_csv ES01" \
        -i "pandas.read_html ES01" \
        -i "pandas.DataFrame.to_html ES01" \
        -i "pandas.io.formats.style.Styler.to_html ES01" \
        -i "pandas.read_xml ES01" \
        -i "pandas.DataFrame.to_xml ES01" \
        -i "pandas.io.formats.style.Styler.to_latex ES01" \
        -i "pandas.HDFStore.get ES01" \
        -i "pandas.HDFStore.info ES01" \
        -i "pandas.HDFStore.keys ES01" \
        -i "pandas.DataFrame.to_feather ES01" \
        -i "pandas.DataFrame.to_orc ES01" \
        -i "pandas.read_sas ES01" \
        -i "pandas.read_spss ES01" \
        -i "pandas.read_stata ES01" \
        -i "pandas.plotting.parallel_coordinates ES01" \
        -i "pandas.bdate_range ES01" \
        -i "pandas.timedelta_range ES01" \
        -i "pandas.interval_range ES01" \
        -i "pandas.util.hash_array ES01" \
        -i "pandas.Index.is_monotonic_increasing ES01" \
        -i "pandas.Index.is_monotonic_decreasing ES01" \
        -i "pandas.Index.is_unique ES01" \
        -i "pandas.Index.has_duplicates ES01" \
        -i "pandas.Index.dtype ES01" \
        -i "pandas.Index.inferred_type ES01" \
        -i "pandas.Index.shape ES01" \
        -i "pandas.Index.name ES01" \
        -i "pandas.Index.memory_usage ES01" \
        -i "pandas.Index.all ES01" \
        -i "pandas.Index.any ES01" \
        -i "pandas.Index.delete ES01" \
        -i "pandas.Index.drop ES01" \
        -i "pandas.Index.drop_duplicates ES01" \
        -i "pandas.Index.identical ES01" \
        -i "pandas.Index.min ES01" \
        -i "pandas.Index.max ES01" \
        -i "pandas.Index.reindex ES01" \
        -i "pandas.Index.putmask ES01" \
        -i "pandas.Index.fillna ES01" \
        -i "pandas.Index.dropna ES01" \
        -i "pandas.Index.infer_objects ES01" \
        -i "pandas.Index.map ES01" \
        -i "pandas.Index.ravel ES01" \
        -i "pandas.Index.argsort ES01" \
        -i "pandas.Index.append ES01" \
        -i "pandas.Index.join ES01" \
        -i "pandas.Index.symmetric_difference ES01" \
        -i "pandas.Index.get_loc ES01" \
        -i "pandas.Index.slice_locs ES01" \
        -i "pandas.CategoricalIndex.append ES01" \
        -i "pandas.IntervalIndex ES01" \
        -i "pandas.IntervalIndex.from_arrays ES01" \
        -i "pandas.IntervalIndex.from_tuples ES01" \
        -i "pandas.IntervalIndex.from_breaks ES01" \
        -i "pandas.IndexSlice ES01" \
        -i "pandas.TimedeltaIndex.to_pytimedelta ES01" \
        -i "pandas.PeriodIndex.day ES01" \
        -i "pandas.PeriodIndex.dayofweek ES01" \
        -i "pandas.PeriodIndex.day_of_week ES01" \
        -i "pandas.PeriodIndex.dayofyear ES01" \
        -i "pandas.PeriodIndex.day_of_year ES01" \
        -i "pandas.PeriodIndex.days_in_month ES01" \
        -i "pandas.PeriodIndex.daysinmonth ES01" \
        -i "pandas.PeriodIndex.hour ES01" \
        -i "pandas.PeriodIndex.is_leap_year ES01" \
        -i "pandas.PeriodIndex.minute ES01" \
        -i "pandas.PeriodIndex.month ES01" \
        -i "pandas.PeriodIndex.quarter ES01" \
        -i "pandas.PeriodIndex.second ES01" \
        -i "pandas.PeriodIndex.week ES01" \
        -i "pandas.PeriodIndex.weekday ES01" \
        -i "pandas.PeriodIndex.weekofyear ES01" \
        -i "pandas.PeriodIndex.year ES01" \
        -i "pandas.PeriodIndex.from_fields ES01" \
        -i "pandas.PeriodIndex.from_ordinals ES01" \
        -i "pandas.api.typing.Rolling.count ES01" \
        -i "pandas.api.typing.Rolling.sum ES01" \
        -i "pandas.api.typing.Rolling.mean ES01" \
        -i "pandas.api.typing.Rolling.median ES01" \
        -i "pandas.api.typing.Rolling.var ES01" \
        -i "pandas.api.typing.Rolling.std ES01" \
        -i "pandas.api.typing.Rolling.min ES01" \
        -i "pandas.api.typing.Rolling.max ES01" \
        -i "pandas.api.typing.Rolling.first ES01" \
        -i "pandas.api.typing.Rolling.last ES01" \
        -i "pandas.api.typing.Rolling.corr ES01" \
        -i "pandas.api.typing.Rolling.cov ES01" \
        -i "pandas.api.typing.Rolling.skew ES01" \
        -i "pandas.api.typing.Rolling.kurt ES01" \
        -i "pandas.api.typing.Rolling.apply ES01" \
        -i "pandas.api.typing.Rolling.aggregate ES01" \
        -i "pandas.api.typing.Rolling.quantile ES01" \
        -i "pandas.api.typing.Rolling.sem ES01" \
        -i "pandas.api.typing.Rolling.rank ES01" \
        -i "pandas.api.typing.Window.mean ES01" \
        -i "pandas.api.typing.Window.sum ES01" \
        -i "pandas.api.typing.Window.var ES01" \
        -i "pandas.api.typing.Window.std ES01" \
        -i "pandas.api.typing.Expanding.count ES01" \
        -i "pandas.api.typing.Expanding.sum ES01" \
        -i "pandas.api.typing.Expanding.mean ES01" \
        -i "pandas.api.typing.Expanding.median ES01" \
        -i "pandas.api.typing.Expanding.var ES01" \
        -i "pandas.api.typing.Expanding.std ES01" \
        -i "pandas.api.typing.Expanding.min ES01" \
        -i "pandas.api.typing.Expanding.max ES01" \
        -i "pandas.api.typing.Expanding.first ES01" \
        -i "pandas.api.typing.Expanding.last ES01" \
        -i "pandas.api.typing.Expanding.corr ES01" \
        -i "pandas.api.typing.Expanding.cov ES01" \
        -i "pandas.api.typing.Expanding.skew ES01" \
        -i "pandas.api.typing.Expanding.kurt ES01" \
        -i "pandas.api.typing.Expanding.apply ES01" \
        -i "pandas.api.typing.Expanding.aggregate ES01" \
        -i "pandas.api.typing.Expanding.quantile ES01" \
        -i "pandas.api.typing.Expanding.sem ES01" \
        -i "pandas.api.typing.Expanding.rank ES01" \
        -i "pandas.api.typing.ExponentialMovingWindow.mean ES01" \
        -i "pandas.api.typing.ExponentialMovingWindow.sum ES01" \
        -i "pandas.api.typing.ExponentialMovingWindow.std ES01" \
        -i "pandas.api.typing.ExponentialMovingWindow.var ES01" \
        -i "pandas.api.typing.ExponentialMovingWindow.corr ES01" \
        -i "pandas.api.typing.ExponentialMovingWindow.cov ES01" \
        -i "pandas.api.indexers.BaseIndexer ES01" \
        -i "pandas.api.indexers.FixedForwardWindowIndexer ES01" \
        -i "pandas.api.indexers.VariableOffsetWindowIndexer ES01" \
        -i "pandas.NamedAgg ES01" \
        -i "pandas.api.typing.DataFrameGroupBy.all ES01" \
        -i "pandas.api.typing.DataFrameGroupBy.any ES01" \
        -i "pandas.api.typing.DataFrameGroupBy.bfill ES01" \
        -i "pandas.api.typing.DataFrameGroupBy.corr ES01" \
        -i "pandas.api.typing.DataFrameGroupBy.count ES01" \
        -i "pandas.api.typing.DataFrameGroupBy.cummin ES01" \
        -i "pandas.api.typing.DataFrameGroupBy.cumprod ES01" \
        -i "pandas.api.typing.DataFrameGroupBy.cumsum ES01" \
        -i "pandas.api.typing.DataFrameGroupBy.ewm ES01" \
        -i "pandas.api.typing.DataFrameGroupBy.expanding ES01" \
        -i "pandas.api.typing.DataFrameGroupBy.ffill ES01" \
        -i "pandas.api.typing.DataFrameGroupBy.idxmax ES01" \
        -i "pandas.api.typing.DataFrameGroupBy.idxmin ES01" \
        -i "pandas.api.typing.DataFrameGroupBy.max ES01" \
        -i "pandas.api.typing.DataFrameGroupBy.mean ES01" \
        -i "pandas.api.typing.DataFrameGroupBy.min ES01" \
        -i "pandas.api.typing.DataFrameGroupBy.nunique ES01" \
        -i "pandas.api.typing.DataFrameGroupBy.pct_change ES01" \
        -i "pandas.api.typing.DataFrameGroupBy.prod ES01" \
        -i "pandas.api.typing.DataFrameGroupBy.quantile ES01" \
        -i "pandas.api.typing.DataFrameGroupBy.rank ES01" \
        -i "pandas.api.typing.DataFrameGroupBy.rolling ES01" \
        -i "pandas.api.typing.DataFrameGroupBy.size ES01" \
        -i "pandas.api.typing.DataFrameGroupBy.kurt ES01" \
        -i "pandas.api.typing.DataFrameGroupBy.sum ES01" \
        -i "pandas.api.typing.DataFrameGroupBy.value_counts ES01" \
        -i "pandas.api.typing.SeriesGroupBy.all ES01" \
        -i "pandas.api.typing.SeriesGroupBy.any ES01" \
        -i "pandas.api.typing.SeriesGroupBy.bfill ES01" \
        -i "pandas.api.typing.SeriesGroupBy.corr ES01" \
        -i "pandas.api.typing.SeriesGroupBy.count ES01" \
        -i "pandas.api.typing.SeriesGroupBy.cov ES01" \
        -i "pandas.api.typing.SeriesGroupBy.cummin ES01" \
        -i "pandas.api.typing.SeriesGroupBy.cumprod ES01" \
        -i "pandas.api.typing.SeriesGroupBy.cumsum ES01" \
        -i "pandas.api.typing.SeriesGroupBy.ewm ES01" \
        -i "pandas.api.typing.SeriesGroupBy.expanding ES01" \
        -i "pandas.api.typing.SeriesGroupBy.ffill ES01" \
        -i "pandas.api.typing.SeriesGroupBy.is_monotonic_increasing ES01" \
        -i "pandas.api.typing.SeriesGroupBy.is_monotonic_decreasing ES01" \
        -i "pandas.api.typing.SeriesGroupBy.max ES01" \
        -i "pandas.api.typing.SeriesGroupBy.mean ES01" \
        -i "pandas.api.typing.SeriesGroupBy.min ES01" \
        -i "pandas.api.typing.SeriesGroupBy.nlargest ES01" \
        -i "pandas.api.typing.SeriesGroupBy.nsmallest ES01" \
        -i "pandas.api.typing.SeriesGroupBy.nunique ES01" \
        -i "pandas.api.typing.SeriesGroupBy.pct_change ES01" \
        -i "pandas.api.typing.SeriesGroupBy.prod ES01" \
        -i "pandas.api.typing.SeriesGroupBy.quantile ES01" \
        -i "pandas.api.typing.SeriesGroupBy.rank ES01" \
        -i "pandas.api.typing.SeriesGroupBy.rolling ES01" \
        -i "pandas.api.typing.SeriesGroupBy.size ES01" \
        -i "pandas.api.typing.SeriesGroupBy.kurt ES01" \
        -i "pandas.api.typing.SeriesGroupBy.sum ES01" \
        -i "pandas.api.typing.SeriesGroupBy.value_counts ES01" \
        -i "pandas.api.typing.DataFrameGroupBy.boxplot ES01" \
        -i "pandas.api.typing.SeriesGroupBy.hist ES01" \
        -i "pandas.io.formats.style.Styler.set_uuid ES01" \
        -i "pandas.io.formats.style.Styler.pipe ES01" \
        -i "pandas.io.formats.style.Styler.highlight_null ES01" \
        -i "pandas.io.formats.style.Styler.highlight_max ES01" \
        -i "pandas.io.formats.style.Styler.highlight_min ES01" \
        -i "pandas.io.formats.style.Styler.highlight_between ES01" \
        -i "pandas.io.formats.style.Styler.highlight_quantile ES01" \
        -i "pandas.io.formats.style.Styler.bar ES01" \
        -i "pandas.io.formats.style.Styler.to_string ES01" \
        -i "pandas.api.extensions.register_dataframe_accessor ES01" \
        -i "pandas.api.extensions.register_series_accessor ES01" \
        -i "pandas.api.extensions.register_index_accessor ES01" \
        -i "pandas.DataFrame.agg ES01" \
        -i "pandas.DataFrame.aggregate ES01" \
        -i "pandas.DataFrame.transform ES01" \
        -i "pandas.DataFrame.rolling ES01" \
        -i "pandas.DataFrame.corr ES01" \
        -i "pandas.DataFrame.mean ES01" \
        -i "pandas.DataFrame.median ES01" \
        -i "pandas.DataFrame.prod ES01" \
        -i "pandas.DataFrame.product ES01" \
        -i "pandas.DataFrame.quantile ES01" \
        -i "pandas.DataFrame.round ES01" \
        -i "pandas.DataFrame.value_counts ES01" \
        -i "pandas.DataFrame.at_time ES01" \
        -i "pandas.DataFrame.to_markdown ES01" # no backslash in the last line

    RET=$(($RET + $?)) ; echo $MSG "DONE"

fi

### DOCUMENTATION NOTEBOOKS ###
if [[ -z "$CHECK" || "$CHECK" == "notebooks" ]]; then

    MSG='Notebooks' ; echo $MSG
    jupyter nbconvert --execute "$(find doc/source -name '*.ipynb')" --to notebook
    RET=$(($RET + $?)) ; echo $MSG "DONE"

fi

### SINGLE-PAGE DOCS ###
if [[ -z "$CHECK" || "$CHECK" == "single-docs" ]]; then
    python doc/make.py --warnings-are-errors --no-browser --single pandas.Series.value_counts
    python doc/make.py --warnings-are-errors --no-browser --single pandas.Series.str.split
fi

exit $RET
