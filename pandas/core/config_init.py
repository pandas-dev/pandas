import pandas.core.config as cf
from pandas.core.config import is_int, is_bool, is_text, is_float
from pandas.core.format import detect_console_encoding

"""
This module is imported from the pandas package __init__.py file
in order to ensure that the core.config options registered here will
be available as soon as the user loads the package. if register_option
is invoked inside specific modules, they will not be registered until that
module is imported, which may or may not be a problem.

If you need to make sure options are available even before a certain
module is imported, register them here rather then in the module.

"""


###########################################
# options from the "display" namespace

pc_precision_doc = """
: int
    Floating point output precision (number of significant digits). This is
    only a suggestion
"""

pc_colspace_doc = """
: int
    Default space for DataFrame columns, defaults to 12
"""

pc_max_rows_doc = """
: int
    This sets the maximum number of rows pandas should output when printing
    out various output. For example, this value determines whether the repr()
    for a dataframe prints out fully or just an summary repr.
"""

pc_max_cols_doc = """
: int
    max_rows and max_columns are used in __repr__() methods to decide if
    to_string() or info() is used to render an object to a string.
    Either one, or both can be set to 0 (experimental). Pandas will figure
    out how big the terminal is and will not display more rows or/and
    columns that can fit on it.
"""

pc_max_info_cols_doc = """
: int
    max_info_columns is used in DataFrame.info method to decide if
    per column information will be printed.
"""

pc_nb_repr_h_doc = """
: boolean
    When True (default), IPython notebook will use html representation for
    pandas objects (if it is available).
"""

pc_date_dayfirst_doc = """
: boolean
    When True, prints and parses dates with the day first, eg 20/01/2005
"""

pc_date_yearfirst_doc = """
: boolean
    When True, prints and parses dates with the year first, eg 2005/01/20
"""

pc_pprint_nest_depth = """
: int
    Defaults to 3.
    Controls the number of nested levels to process when pretty-printing
"""

pc_multi_sparse_doc = """
: boolean
    Default True, "sparsify" MultiIndex display (don't display repeated
    elements in outer levels within groups)
"""

pc_encoding_doc = """
: str/unicode
    Defaults to the detected encoding of the console.
    Specifies the encoding to be used for strings returned by to_string,
    these are generally strings meant to be displayed on the console.
"""

float_format_doc = """
: callable
    The callable should accept a floating point number and return
    a string with the desired format of the number. This is used
    in some places like SeriesFormatter.
    See core.format.EngFormatter for an example.

"""

max_colwidth_doc = """
: int
    The maximum width in characters of a column in the repr of
    a pandas data structure. When the column overflows, a "..."
    placeholder is embedded in the output.
"""

colheader_justify_doc = """
: 'left'/'right'
    Controls the justification of column headers. used by DataFrameFormatter.
"""

pc_expand_repr_doc = """
: boolean
    Default False
    Whether to print out the full DataFrame repr for wide DataFrames
    across multiple lines.
    If False, the summary representation is shown.
"""

pc_line_width_doc = """
: int
    Default 80
    When printing wide DataFrames, this is the width of each line.
"""

with cf.config_prefix('display'):
    cf.register_option('precision', 7, pc_precision_doc, validator=is_int)
    cf.register_option('float_format', None, float_format_doc)
    cf.register_option('column_space', 12, validator=is_int)
    cf.register_option('max_rows', 100, pc_max_rows_doc, validator=is_int)
    cf.register_option('max_colwidth', 50, max_colwidth_doc, validator=is_int)
    cf.register_option('max_columns', 20, pc_max_cols_doc, validator=is_int)
    cf.register_option('max_info_columns', 100, pc_max_info_cols_doc,
                       validator=is_int)
    cf.register_option('colheader_justify', 'right', colheader_justify_doc,
                       validator=is_text)
    cf.register_option('notebook_repr_html', True, pc_nb_repr_h_doc,
                       validator=is_bool)
    cf.register_option('date_dayfirst', False, pc_date_dayfirst_doc,
                       validator=is_bool)
    cf.register_option('date_yearfirst', False, pc_date_yearfirst_doc,
                       validator=is_bool)
    cf.register_option('pprint_nest_depth', 3, pc_pprint_nest_depth,
                       validator=is_int)
    cf.register_option('multi_sparse', True, pc_multi_sparse_doc,
                       validator=is_bool)
    cf.register_option('encoding', detect_console_encoding(), pc_encoding_doc,
                       validator=is_text)
    cf.register_option('expand_frame_repr', True, pc_expand_repr_doc)
    cf.register_option('line_width', 80, pc_line_width_doc)

tc_sim_interactive_doc = """
: boolean
    Default False
    Whether to simulate interactive mode for purposes of testing
"""
with cf.config_prefix('mode'):
    cf.register_option('sim_interactive', False, tc_sim_interactive_doc)

use_inf_as_null_doc = """
: boolean
    True means treat None, NaN, INF, -INF as null (old way),
    False means None and NaN are null, but INF, -INF are not null
    (new way).
"""

# we don't want to start importing evrything at the global context level
# or we'll hit circular deps.


def use_inf_as_null_cb(key):
    from pandas.core.common import _use_inf_as_null
    _use_inf_as_null(key)

with cf.config_prefix('mode'):
    cf.register_option('use_inf_as_null', False, use_inf_as_null_doc,
                       cb=use_inf_as_null_cb)
