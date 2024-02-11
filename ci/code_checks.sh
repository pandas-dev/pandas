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
    { echo "Unknown command $1. Usage: $0 [code|doctests|docstrings|single-docs|notebooks]"; exit 9999; }

BASE_DIR="$(dirname $0)/.."
RET=0

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

    MSG='Python and Cython Doctests' ; echo $MSG
    python -c 'import pandas as pd; pd.test(run_doctests=True)'
    RET=$(($RET + $?)) ; echo $MSG "DONE"

fi

### DOCSTRINGS ###
if [[ -z "$CHECK" || "$CHECK" == "docstrings" ]]; then

    MSG='Validate docstrings (EX01, EX03, EX04, GL01, GL02, GL03, GL04, GL05, GL06, GL07, GL09, GL10, PR03, PR04, PR05, PR06, PR08, PR09, PR10, RT01, RT02, RT04, RT05, SA02, SA03, SA04, SS01, SS02, SS03, SS04, SS05, SS06)' ; echo $MSG
    $BASE_DIR/scripts/validate_docstrings.py --format=actions --errors=EX01,EX03,EX04,GL01,GL02,GL03,GL04,GL05,GL06,GL07,GL09,GL10,PR03,PR04,PR05,PR06,PR08,PR09,PR10,RT01,RT02,RT04,RT05,SA02,SA03,SA04,SS01,SS02,SS03,SS04,SS05,SS06
    RET=$(($RET + $?)) ; echo $MSG "DONE"

    MSG='Partially validate docstrings (PR02)' ;  echo $MSG
    $BASE_DIR/scripts/validate_docstrings.py --format=actions --errors=PR02 --ignore_functions \
        pandas.Series.dt.to_period\
        pandas.Series.dt.tz_localize\
        pandas.Series.dt.tz_convert\
        pandas.Series.dt.strftime\
        pandas.Series.dt.round\
        pandas.Series.dt.floor\
        pandas.Series.dt.ceil\
        pandas.Series.dt.month_name\
        pandas.Series.dt.day_name\
        pandas.Series.str.wrap\
        pandas.Series.cat.rename_categories\
        pandas.Series.cat.reorder_categories\
        pandas.Series.cat.add_categories\
        pandas.Series.cat.remove_categories\
        pandas.Series.cat.set_categories\
        pandas.Series.plot\
        pandas.DataFrame.plot\
        pandas.tseries.offsets.DateOffset\
        pandas.tseries.offsets.BusinessDay\
        pandas.tseries.offsets.BDay\
        pandas.tseries.offsets.BusinessHour\
        pandas.tseries.offsets.CustomBusinessDay\
        pandas.tseries.offsets.CDay\
        pandas.tseries.offsets.CustomBusinessHour\
        pandas.tseries.offsets.MonthEnd\
        pandas.tseries.offsets.MonthBegin\
        pandas.tseries.offsets.BusinessMonthEnd\
        pandas.tseries.offsets.BMonthEnd\
        pandas.tseries.offsets.BusinessMonthBegin\
        pandas.tseries.offsets.BMonthBegin\
        pandas.tseries.offsets.CustomBusinessMonthEnd\
        pandas.tseries.offsets.CBMonthEnd\
        pandas.tseries.offsets.CustomBusinessMonthBegin\
        pandas.tseries.offsets.CBMonthBegin\
        pandas.tseries.offsets.SemiMonthEnd\
        pandas.tseries.offsets.SemiMonthBegin\
        pandas.tseries.offsets.Week\
        pandas.tseries.offsets.WeekOfMonth\
        pandas.tseries.offsets.LastWeekOfMonth\
        pandas.tseries.offsets.BQuarterEnd\
        pandas.tseries.offsets.BQuarterBegin\
        pandas.tseries.offsets.QuarterEnd\
        pandas.tseries.offsets.QuarterBegin\
        pandas.tseries.offsets.BYearEnd\
        pandas.tseries.offsets.BYearBegin\
        pandas.tseries.offsets.YearEnd\
        pandas.tseries.offsets.YearBegin\
        pandas.tseries.offsets.FY5253\
        pandas.tseries.offsets.FY5253Quarter\
        pandas.tseries.offsets.Easter\
        pandas.tseries.offsets.Day\
        pandas.tseries.offsets.Hour\
        pandas.tseries.offsets.Minute\
        pandas.tseries.offsets.Second\
        pandas.tseries.offsets.Milli\
        pandas.tseries.offsets.Micro\
        pandas.tseries.offsets.Nano\
        pandas.Timestamp.max\
        pandas.Timestamp.min\
        pandas.Timestamp.resolution\
        pandas.Timedelta.max\
        pandas.Timedelta.min\
        pandas.Timedelta.resolution\
        pandas.Interval\
        pandas.Grouper\
        pandas.core.groupby.DataFrameGroupBy.nth\
        pandas.core.groupby.DataFrameGroupBy.rolling\
        pandas.core.groupby.SeriesGroupBy.nth\
        pandas.core.groupby.SeriesGroupBy.rolling\
        pandas.core.groupby.DataFrameGroupBy.plot\
        pandas.core.groupby.SeriesGroupBy.plot # There should be no backslash in the final line, please keep this comment in the last ignored function
    RET=$(($RET + $?)) ; echo $MSG "DONE"

    MSG='Partially validate docstrings (RT03)' ;  echo $MSG
    $BASE_DIR/scripts/validate_docstrings.py --format=actions --errors=RT03 --ignore_functions \
        pandas.Categorical.set_categories\
        pandas.CategoricalIndex.set_categories\
        pandas.DataFrame.astype\
        pandas.DataFrame.at_time\
        pandas.DataFrame.ewm\
        pandas.DataFrame.expanding\
        pandas.DataFrame.filter\
        pandas.DataFrame.first_valid_index\
        pandas.DataFrame.get\
        pandas.DataFrame.hist\
        pandas.DataFrame.infer_objects\
        pandas.DataFrame.kurt\
        pandas.DataFrame.kurtosis\
        pandas.DataFrame.last_valid_index\
        pandas.DataFrame.mask\
        pandas.DataFrame.max\
        pandas.DataFrame.mean\
        pandas.DataFrame.median\
        pandas.DataFrame.min\
        pandas.DataFrame.nsmallest\
        pandas.DataFrame.nunique\
        pandas.DataFrame.pipe\
        pandas.DataFrame.plot.box\
        pandas.DataFrame.plot.density\
        pandas.DataFrame.plot.kde\
        pandas.DataFrame.plot.scatter\
        pandas.DataFrame.pop\
        pandas.DataFrame.prod\
        pandas.DataFrame.product\
        pandas.DataFrame.reindex\
        pandas.DataFrame.reorder_levels\
        pandas.DataFrame.sem\
        pandas.DataFrame.skew\
        pandas.DataFrame.std\
        pandas.DataFrame.sum\
        pandas.DataFrame.swapaxes\
        pandas.DataFrame.to_numpy\
        pandas.DataFrame.to_orc\
        pandas.DataFrame.to_parquet\
        pandas.DataFrame.unstack\
        pandas.DataFrame.value_counts\
        pandas.DataFrame.var\
        pandas.DataFrame.where\
        pandas.DatetimeIndex.indexer_at_time\
        pandas.DatetimeIndex.indexer_between_time\
        pandas.DatetimeIndex.snap\
        pandas.DatetimeIndex.std\
        pandas.DatetimeIndex.to_period\
        pandas.DatetimeIndex.to_pydatetime\
        pandas.DatetimeIndex.tz_convert\
        pandas.HDFStore.info\
        pandas.Index.append\
        pandas.Index.difference\
        pandas.Index.drop_duplicates\
        pandas.Index.droplevel\
        pandas.Index.dropna\
        pandas.Index.duplicated\
        pandas.Index.fillna\
        pandas.Index.get_loc\
        pandas.Index.insert\
        pandas.Index.intersection\
        pandas.Index.join\
        pandas.Index.memory_usage\
        pandas.Index.nunique\
        pandas.Index.putmask\
        pandas.Index.ravel\
        pandas.Index.slice_indexer\
        pandas.Index.slice_locs\
        pandas.Index.symmetric_difference\
        pandas.Index.to_list\
        pandas.Index.union\
        pandas.Index.unique\
        pandas.Index.value_counts\
        pandas.IntervalIndex.contains\
        pandas.IntervalIndex.get_loc\
        pandas.IntervalIndex.set_closed\
        pandas.IntervalIndex.to_tuples\
        pandas.MultiIndex.copy\
        pandas.MultiIndex.drop\
        pandas.MultiIndex.droplevel\
        pandas.MultiIndex.remove_unused_levels\
        pandas.MultiIndex.reorder_levels\
        pandas.MultiIndex.set_levels\
        pandas.MultiIndex.to_frame\
        pandas.PeriodIndex.to_timestamp\
        pandas.Series.__iter__\
        pandas.Series.astype\
        pandas.Series.at_time\
        pandas.Series.case_when\
        pandas.Series.cat.set_categories\
        pandas.Series.dt.to_period\
        pandas.Series.dt.tz_convert\
        pandas.Series.ewm\
        pandas.Series.expanding\
        pandas.Series.filter\
        pandas.Series.first_valid_index\
        pandas.Series.get\
        pandas.Series.infer_objects\
        pandas.Series.kurt\
        pandas.Series.kurtosis\
        pandas.Series.last_valid_index\
        pandas.Series.mask\
        pandas.Series.max\
        pandas.Series.mean\
        pandas.Series.median\
        pandas.Series.min\
        pandas.Series.nunique\
        pandas.Series.pipe\
        pandas.Series.plot.box\
        pandas.Series.plot.density\
        pandas.Series.plot.kde\
        pandas.Series.pop\
        pandas.Series.prod\
        pandas.Series.product\
        pandas.Series.reindex\
        pandas.Series.reorder_levels\
        pandas.Series.sem\
        pandas.Series.skew\
        pandas.Series.sparse.to_coo\
        pandas.Series.std\
        pandas.Series.str.capitalize\
        pandas.Series.str.casefold\
        pandas.Series.str.center\
        pandas.Series.str.decode\
        pandas.Series.str.encode\
        pandas.Series.str.find\
        pandas.Series.str.fullmatch\
        pandas.Series.str.get\
        pandas.Series.str.index\
        pandas.Series.str.ljust\
        pandas.Series.str.lower\
        pandas.Series.str.lstrip\
        pandas.Series.str.match\
        pandas.Series.str.normalize\
        pandas.Series.str.partition\
        pandas.Series.str.rfind\
        pandas.Series.str.rindex\
        pandas.Series.str.rjust\
        pandas.Series.str.rpartition\
        pandas.Series.str.rstrip\
        pandas.Series.str.strip\
        pandas.Series.str.swapcase\
        pandas.Series.str.title\
        pandas.Series.str.translate\
        pandas.Series.str.upper\
        pandas.Series.str.wrap\
        pandas.Series.str.zfill\
        pandas.Series.sum\
        pandas.Series.to_list\
        pandas.Series.to_numpy\
        pandas.Series.to_timestamp\
        pandas.Series.value_counts\
        pandas.Series.var\
        pandas.Series.where\
        pandas.TimedeltaIndex.as_unit\
        pandas.TimedeltaIndex.to_pytimedelta\
        pandas.api.extensions.ExtensionArray._accumulate\
        pandas.api.extensions.ExtensionArray._hash_pandas_object\
        pandas.api.extensions.ExtensionArray._pad_or_backfill\
        pandas.api.extensions.ExtensionArray._reduce\
        pandas.api.extensions.ExtensionArray.copy\
        pandas.api.extensions.ExtensionArray.dropna\
        pandas.api.extensions.ExtensionArray.duplicated\
        pandas.api.extensions.ExtensionArray.insert\
        pandas.api.extensions.ExtensionArray.isin\
        pandas.api.extensions.ExtensionArray.ravel\
        pandas.api.extensions.ExtensionArray.take\
        pandas.api.extensions.ExtensionArray.tolist\
        pandas.api.extensions.ExtensionArray.unique\
        pandas.api.interchange.from_dataframe\
        pandas.api.types.is_hashable\
        pandas.api.types.pandas_dtype\
        pandas.api.types.union_categoricals\
        pandas.arrays.IntervalArray.contains\
        pandas.arrays.IntervalArray.set_closed\
        pandas.arrays.IntervalArray.to_tuples\
        pandas.bdate_range\
        pandas.core.groupby.DataFrameGroupBy.__iter__\
        pandas.core.groupby.DataFrameGroupBy.agg\
        pandas.core.groupby.DataFrameGroupBy.aggregate\
        pandas.core.groupby.DataFrameGroupBy.apply\
        pandas.core.groupby.DataFrameGroupBy.boxplot\
        pandas.core.groupby.DataFrameGroupBy.cummax\
        pandas.core.groupby.DataFrameGroupBy.cummin\
        pandas.core.groupby.DataFrameGroupBy.cumprod\
        pandas.core.groupby.DataFrameGroupBy.cumsum\
        pandas.core.groupby.DataFrameGroupBy.filter\
        pandas.core.groupby.DataFrameGroupBy.get_group\
        pandas.core.groupby.DataFrameGroupBy.hist\
        pandas.core.groupby.DataFrameGroupBy.mean\
        pandas.core.groupby.DataFrameGroupBy.nunique\
        pandas.core.groupby.DataFrameGroupBy.rank\
        pandas.core.groupby.DataFrameGroupBy.resample\
        pandas.core.groupby.DataFrameGroupBy.skew\
        pandas.core.groupby.DataFrameGroupBy.transform\
        pandas.core.groupby.SeriesGroupBy.__iter__\
        pandas.core.groupby.SeriesGroupBy.agg\
        pandas.core.groupby.SeriesGroupBy.aggregate\
        pandas.core.groupby.SeriesGroupBy.apply\
        pandas.core.groupby.SeriesGroupBy.cummax\
        pandas.core.groupby.SeriesGroupBy.cummin\
        pandas.core.groupby.SeriesGroupBy.cumprod\
        pandas.core.groupby.SeriesGroupBy.cumsum\
        pandas.core.groupby.SeriesGroupBy.filter\
        pandas.core.groupby.SeriesGroupBy.get_group\
        pandas.core.groupby.SeriesGroupBy.mean\
        pandas.core.groupby.SeriesGroupBy.rank\
        pandas.core.groupby.SeriesGroupBy.resample\
        pandas.core.groupby.SeriesGroupBy.skew\
        pandas.core.groupby.SeriesGroupBy.transform\
        pandas.core.resample.Resampler.__iter__\
        pandas.core.resample.Resampler.ffill\
        pandas.core.resample.Resampler.get_group\
        pandas.core.resample.Resampler.max\
        pandas.core.resample.Resampler.min\
        pandas.core.resample.Resampler.transform\
        pandas.date_range\
        pandas.interval_range\
        pandas.io.formats.style.Styler.apply\
        pandas.io.formats.style.Styler.apply_index\
        pandas.io.formats.style.Styler.background_gradient\
        pandas.io.formats.style.Styler.bar\
        pandas.io.formats.style.Styler.concat\
        pandas.io.formats.style.Styler.export\
        pandas.io.formats.style.Styler.format\
        pandas.io.formats.style.Styler.format_index\
        pandas.io.formats.style.Styler.hide\
        pandas.io.formats.style.Styler.highlight_between\
        pandas.io.formats.style.Styler.highlight_max\
        pandas.io.formats.style.Styler.highlight_min\
        pandas.io.formats.style.Styler.highlight_null\
        pandas.io.formats.style.Styler.highlight_quantile\
        pandas.io.formats.style.Styler.map\
        pandas.io.formats.style.Styler.map_index\
        pandas.io.formats.style.Styler.relabel_index\
        pandas.io.formats.style.Styler.set_caption\
        pandas.io.formats.style.Styler.set_properties\
        pandas.io.formats.style.Styler.set_sticky\
        pandas.io.formats.style.Styler.set_table_attributes\
        pandas.io.formats.style.Styler.set_table_styles\
        pandas.io.formats.style.Styler.set_td_classes\
        pandas.io.formats.style.Styler.set_tooltips\
        pandas.io.formats.style.Styler.set_uuid\
        pandas.io.formats.style.Styler.text_gradient\
        pandas.io.formats.style.Styler.use\
        pandas.io.json.build_table_schema\
        pandas.io.stata.StataReader.value_labels\
        pandas.io.stata.StataReader.variable_labels\
        pandas.json_normalize\
        pandas.merge_asof\
        pandas.period_range\
        pandas.plotting.andrews_curves\
        pandas.plotting.autocorrelation_plot\
        pandas.plotting.lag_plot\
        pandas.plotting.parallel_coordinates\
        pandas.plotting.radviz\
        pandas.plotting.table\
        pandas.read_feather\
        pandas.read_orc\
        pandas.read_sas\
        pandas.read_spss\
        pandas.read_sql\
        pandas.read_sql_query\
        pandas.read_stata\
        pandas.set_eng_float_format\
        pandas.timedelta_range\
        pandas.util.hash_pandas_object # There should be no backslash in the final line, please keep this comment in the last ignored function
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
    python doc/make.py --warnings-are-errors --no-browser --single pandas.Series.value_counts
    python doc/make.py --warnings-are-errors --no-browser --single pandas.Series.str.split
fi

exit $RET
