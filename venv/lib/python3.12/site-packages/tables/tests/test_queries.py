"""Test module for queries on datasets."""

import re
import sys
import warnings
import functools

import numpy as np

import tables as tb
from tables.tests import common

# Data parameters
# ---------------
row_period = 50
"""Maximum number of unique rows before they start cycling."""
md_shape = (2, 2)
"""Shape of multidimensional fields."""

_maxnvalue = row_period + np.prod(md_shape, dtype=tb.utils.SizeType) - 1
_strlen = int(np.log10(_maxnvalue - 1)) + 1

str_format = "%%0%dd" % _strlen
"""Format of string values."""

small_blocksizes = (300, 60, 20, 5)
# small_blocksizes = (512, 128, 32, 4)   # for manual testing only
"""Sensible parameters for indexing with small blocksizes."""


# Type information
# ----------------
type_info = {
    "bool": (np.bool_, bool),
    "int8": (np.int8, int),
    "uint8": (np.uint8, int),
    "int16": (np.int16, int),
    "uint16": (np.uint16, int),
    "int32": (np.int32, int),
    "uint32": (np.uint32, int),
    "int64": (np.int64, int),
    "uint64": (np.uint64, int),
    "float32": (np.float32, float),
    "float64": (np.float64, float),
    "complex64": (np.complex64, complex),
    "complex128": (np.complex128, complex),
    "time32": (np.int32, int),
    "time64": (np.float64, float),
    "enum": (np.uint8, int),  # just for these tests
    "string": ("S%s" % _strlen, np.bytes_),  # just for these tests
}
"""NumPy and Numexpr type for each PyTables type that will be tested."""

# globals dict for eval()
func_info = {
    "log10": np.log10,
    "log": np.log,
    "exp": np.exp,
    "abs": np.abs,
    "sqrt": np.sqrt,
    "sin": np.sin,
    "cos": np.cos,
    "tan": np.tan,
    "arcsin": np.arcsin,
    "arccos": np.arccos,
    "arctan": np.arctan,
}
"""functions and NumPy.ufunc() for each function that will be tested."""


if hasattr(np, "float16"):
    type_info["float16"] = (np.float16, float)
# if hasattr(numpy, 'float96'):
#    type_info['float96'] = (np.float96, float)
# if hasattr(numpy, 'float128'):
#    type_info['float128'] = (np.float128, float)
# if hasattr(numpy, 'complex192'):
#    type_info['complex192'] = (np.complex192, complex)
# if hasattr(numpy, 'complex256'):
#    type_info['complex256'] = (np.complex256, complex)

sctype_from_type = {type_: info[0] for (type_, info) in type_info.items()}
"""Maps PyTables types to NumPy scalar types."""
nxtype_from_type = {type_: info[1] for (type_, info) in type_info.items()}
"""Maps PyTables types to Numexpr types."""

heavy_types = frozenset(["uint8", "int16", "uint16", "float32", "complex64"])
"""PyTables types to be tested only in heavy mode."""

enum = tb.Enum({"n%d" % i: i for i in range(_maxnvalue)})
"""Enumerated type to be used in tests."""


# Table description
# -----------------
def append_columns(classdict, shape=()):
    """Append a ``Col`` of each PyTables data type to the `classdict`.

    A column of a certain TYPE gets called ``c_TYPE``.  The number of
    added columns is returned.

    """
    heavy = common.heavy
    for itype, type_ in enumerate(sorted(type_info)):
        if not heavy and type_ in heavy_types:
            continue  # skip heavy type in non-heavy mode
        colpos = itype + 1
        colname = "c_%s" % type_
        if type_ == "enum":
            base = tb.Atom.from_sctype(sctype_from_type[type_])
            col = tb.EnumCol(enum, enum(0), base, shape=shape, pos=colpos)
        else:
            sctype = sctype_from_type[type_]
            dtype = np.dtype((sctype, shape))
            col = tb.Col.from_dtype(dtype, pos=colpos)
        classdict[colname] = col
    ncols = colpos
    return ncols


def nested_description(classname, pos, shape=()):
    """Return a nested column description with all PyTables data types.

    A column of a certain TYPE gets called ``c_TYPE``.  The nested
    column will be placed in the position indicated by `pos`.

    """
    classdict = {}
    append_columns(classdict, shape=shape)
    classdict["_v_pos"] = pos
    return type(classname, (tb.IsDescription,), classdict)


def table_description(classname, nclassname, shape=()):
    """Return a table description for testing queries.

    The description consists of all PyTables data types, both in the
    top level and in the ``c_nested`` nested column.  A column of a
    certain TYPE gets called ``c_TYPE``.  An extra integer column
    ``c_extra`` is also provided.  If a `shape` is given, it will be
    used for all columns.  Finally, an extra indexed column
    ``c_idxextra`` is added as well in order to provide some basic
    tests for multi-index queries.

    """
    classdict = {}
    colpos = append_columns(classdict, shape)

    ndescr = nested_description(nclassname, colpos, shape=shape)
    classdict["c_nested"] = ndescr
    colpos += 1

    extracol = tb.IntCol(shape=shape, pos=colpos)
    classdict["c_extra"] = extracol
    colpos += 1

    idxextracol = tb.IntCol(shape=shape, pos=colpos)
    classdict["c_idxextra"] = idxextracol
    colpos += 1

    return type(classname, (tb.IsDescription,), classdict)


TableDescription = table_description("TableDescription", "NestedDescription")
"""Unidimensional table description for testing queries."""

MDTableDescription = table_description(
    "MDTableDescription", "MDNestedDescription", shape=md_shape
)
"""Multidimensional table description for testing queries."""


# Table data
# ----------
table_data = {}
"""Cached table data for a given shape and number of rows."""
# Data is cached because computing it row by row is quite slow.  Hop!


def fill_table(table, shape, nrows):
    """Fill the given `table` with `nrows` rows of data.

    Values in the i-th row (where 0 <= i < `row_period`) for a
    multidimensional field with M elements span from i to i + M-1.  For
    subsequent rows, values repeat cyclically.

    The same goes for the ``c_extra`` column, but values range from
    -`row_period`/2 to +`row_period`/2.

    """
    # Reuse already computed data if possible.
    tdata = table_data.get((shape, nrows))
    if tdata is not None:
        table.append(tdata)
        table.flush()
        return

    heavy = common.heavy
    size = int(np.prod(shape, dtype=tb.utils.SizeType))

    row, value = table.row, 0
    for nrow in range(nrows):
        data = np.arange(value, value + size).reshape(shape)
        for type_, sctype in sctype_from_type.items():
            if not heavy and type_ in heavy_types:
                continue  # skip heavy type in non-heavy mode
            colname = "c_%s" % type_
            ncolname = "c_nested/%s" % colname
            if type_ == "bool":
                coldata = data > (row_period // 2)
            elif type_ == "string":
                sdata = [str_format % x for x in range(value, value + size)]
                coldata = np.array(sdata, dtype=sctype).reshape(shape)
            else:
                coldata = np.asarray(data, dtype=sctype)
            row[ncolname] = row[colname] = coldata
            row["c_extra"] = data - (row_period // 2)
            row["c_idxextra"] = data - (row_period // 2)
        row.append()
        value += 1
        if value == row_period:
            value = 0
    table.flush()

    # Make computed data reusable.
    tdata = table.read()
    table_data[(shape, nrows)] = tdata


class SilentlySkipTest(common.unittest.SkipTest):
    pass


# Base test cases
# ---------------
class BaseTableQueryTestCase(common.TempFileMixin, common.PyTablesTestCase):
    """Base test case for querying tables.

    Sub-classes must define the following attributes:

    ``tableDescription``
        The description of the table to be created.
    ``shape``
        The shape of data fields in the table.
    ``nrows``
        The number of data rows to be generated for the table.

    Sub-classes may redefine the following attributes:

    ``indexed``
        Whether columns shall be indexed, if possible.  Default is not
        to index them.
    ``optlevel``
        The level of optimisation of column indexes.  Default is 0.

    """

    indexed = False
    optlevel = 0

    colNotIndexable_re = re.compile(r"\bcan not be indexed\b")
    condNotBoolean_re = re.compile(r"\bdoes not have a boolean type\b")

    def create_indexes(self, colname, ncolname, extracolname):
        if not self.indexed:
            return
        try:
            kind = self.kind
            common.verbosePrint(
                f"* Indexing ``{colname}`` columns. Type: {kind}."
            )
            for acolname in [colname, ncolname, extracolname]:
                acolumn = self.table.colinstances[acolname]
                acolumn.create_index(
                    kind=self.kind,
                    optlevel=self.optlevel,
                    _blocksizes=small_blocksizes,
                    _testmode=True,
                )

        except TypeError as te:
            if self.colNotIndexable_re.search(str(te)):
                raise SilentlySkipTest(
                    "Columns of this type can not be indexed."
                )
            raise
        except NotImplementedError:
            raise SilentlySkipTest(
                "Indexing columns of this type is not supported yet."
            )

    def setUp(self):
        super().setUp()
        self.table = self.h5file.create_table(
            "/", "test", self.tableDescription, expectedrows=self.nrows
        )
        fill_table(self.table, self.shape, self.nrows)


class ScalarTableMixin:
    tableDescription = TableDescription
    shape = ()


class MDTableMixin:
    tableDescription = MDTableDescription
    shape = md_shape


# Test cases on query data
# ------------------------
operators = [
    None,
    "~",
    "<",
    "<=",
    "==",
    "!=",
    ">=",
    ">",
    ("<", "<="),
    (">", ">="),
]
"""Comparison operators to check with different types."""
heavy_operators = frozenset(["~", "<=", ">=", ">", (">", ">=")])
"""Comparison operators to be tested only in heavy mode."""
left_bound = row_period // 4
"""Operand of left side operator in comparisons with operator pairs."""
right_bound = row_period * 3 // 4
"""Operand of right side operator in comparisons with operator pairs."""
func_bound = 0.8  # must be <1 for trig functions to be able to fail
"""Operand of right side operator in comparisons with functions. """
extra_conditions = [
    "",  # uses one index
    "& ((c_extra + 1) < 0)",  # uses one index
    "| (c_idxextra > 0)",  # uses two indexes
    "| ((c_idxextra > 0) | ((c_extra + 1) > 0))",  # can't use indexes
]
"""Extra conditions to append to comparison conditions."""


class TableDataTestCase(BaseTableQueryTestCase):
    """Base test case for querying table data.

    Automatically created test method names have the format
    ``test_XNNNN``, where ``NNNN`` is the zero-padded test number and
    ``X`` indicates whether the test belongs to the light (``l``) or
    heavy (``h``) set.

    """

    _testfmt_light = "test_l%04d"
    _testfmt_heavy = "test_h%04d"


def _old_repr(o):
    if isinstance(o, np.bytes_):
        return repr(bytes(o))
    return repr(o)


def create_test_method(type_, op, extracond, func=None):
    sctype = sctype_from_type[type_]

    # Compute the value of bounds.
    condvars = {
        "bound": right_bound,
        "lbound": left_bound,
        "rbound": right_bound,
        "func_bound": func_bound,
    }
    for bname, bvalue in condvars.items():
        if type_ == "string":
            bvalue = str_format % bvalue
        bvalue = nxtype_from_type[type_](bvalue)
        condvars[bname] = bvalue

    # Compute the name of columns.
    colname = "c_%s" % type_
    ncolname = "c_nested/%s" % colname

    # Compute the query condition.
    if not op:  # as is
        cond = colname
    elif op == "~":  # unary
        cond = "~(%s)" % colname
    elif op == "<" and func is None:  # binary variable-constant
        cond = f'{colname} {op} {_old_repr(condvars["bound"])}'
    elif isinstance(op, tuple):  # double binary variable-constant
        cond = f"(lbound {op[0]} {colname}) & ({colname} {op[1]} rbound)"
    elif func is not None:
        cond = f"{func}({colname}) {op} func_bound"
    else:  # function or binary variable-variable
        cond = f"{colname} {op} bound"
    if extracond:
        cond = f"({cond}) {extracond}"

    def ignore_skipped(oldmethod):
        @functools.wraps(oldmethod)
        def newmethod(self, *args, **kwargs):
            self._verboseHeader()
            try:
                return oldmethod(self, *args, **kwargs)
            except SilentlySkipTest as se:
                if se.args:
                    msg = se.args[0]
                else:
                    msg = "<skipped>"
                common.verbosePrint("\nSkipped test: %s" % msg)
            finally:
                common.verbosePrint("")  # separator line between tests

        return newmethod

    @ignore_skipped
    def test_method(self):
        common.verbosePrint("* Condition is ``%s``." % cond)
        # Replace bitwise operators with their logical counterparts.
        pycond = cond
        for ptop, pyop in [("&", "and"), ("|", "or"), ("~", "not")]:
            pycond = pycond.replace(ptop, pyop)
        pycond = compile(pycond, "<string>", "eval")

        table = self.table
        self.create_indexes(colname, ncolname, "c_idxextra")

        table_slice = dict(start=1, stop=table.nrows - 5, step=3)
        rownos, fvalues = None, None
        # Test that both simple and nested columns work as expected.
        # Knowing how the table is filled, results must be the same.
        for acolname in [colname, ncolname]:
            # First the reference Python version.
            pyrownos, pyfvalues, pyvars = [], [], condvars.copy()
            for row in table.iterrows(**table_slice):
                pyvars[colname] = row[acolname]
                pyvars["c_extra"] = row["c_extra"]
                pyvars["c_idxextra"] = row["c_idxextra"]
                try:
                    with warnings.catch_warnings():
                        warnings.filterwarnings(
                            "ignore",
                            "invalid value encountered in arc(cos|sin)",
                            RuntimeWarning,
                        )
                        isvalidrow = eval(pycond, func_info, pyvars)
                except TypeError:
                    raise SilentlySkipTest(
                        "The Python type does not support the operation."
                    )
                if isvalidrow:
                    pyrownos.append(row.nrow)
                    pyfvalues.append(row[acolname])
            pyrownos = np.array(pyrownos)  # row numbers already sorted
            pyfvalues = np.array(pyfvalues, dtype=sctype)
            pyfvalues.sort()
            common.verbosePrint(
                f"* {len(pyrownos)} rows selected by Python "
                f"from ``{acolname}``."
            )
            if rownos is None:
                rownos = pyrownos  # initialise reference results
                fvalues = pyfvalues
            else:
                self.assertTrue(np.all(pyrownos == rownos))  # check
                self.assertTrue(np.all(pyfvalues == fvalues))

            # Then the in-kernel or indexed version.
            ptvars = condvars.copy()
            ptvars[colname] = table.colinstances[acolname]
            ptvars["c_extra"] = table.colinstances["c_extra"]
            ptvars["c_idxextra"] = table.colinstances["c_idxextra"]
            try:
                isidxq = table.will_query_use_indexing(cond, ptvars)
                # Query twice to trigger possible query result caching.
                ptrownos = [
                    table.get_where_list(
                        cond, condvars, sort=True, **table_slice
                    )
                    for _ in range(2)
                ]
                ptfvalues = [
                    table.read_where(
                        cond, condvars, field=acolname, **table_slice
                    )
                    for _ in range(2)
                ]
            except TypeError as te:
                if self.condNotBoolean_re.search(str(te)):
                    raise SilentlySkipTest("The condition is not boolean.")
                raise
            except NotImplementedError:
                raise SilentlySkipTest(
                    "The PyTables type does not support the operation."
                )
            for ptfvals in ptfvalues:  # row numbers already sorted
                ptfvals.sort()
            common.verbosePrint(
                f"* {len(ptrownos[0])} rows selected by "
                f"PyTables from ``{acolname}``",
                nonl=True,
            )
            common.verbosePrint(f"(indexing: {'yes' if isidxq else 'no'}).")
            self.assertTrue(np.all(ptrownos[0] == rownos))
            self.assertTrue(np.all(ptfvalues[0] == fvalues))
            # The following test possible caching of query results.
            self.assertTrue(np.all(ptrownos[0] == ptrownos[1]))
            self.assertTrue(np.all(ptfvalues[0] == ptfvalues[1]))

    test_method.__doc__ = "Testing ``%s``." % cond
    return test_method


def add_test_method(type_, op, extracond="", func=None):
    global testn
    # Decide to which set the test belongs.
    heavy = type_ in heavy_types or op in heavy_operators
    if heavy:
        testfmt = TableDataTestCase._testfmt_heavy
    else:
        testfmt = TableDataTestCase._testfmt_light
    tmethod = create_test_method(type_, op, extracond, func)
    # The test number is appended to the docstring to help
    # identify failing methods in non-verbose mode.
    tmethod.__name__ = testfmt % testn
    tmethod.__doc__ += testfmt % testn
    setattr(TableDataTestCase, tmethod.__name__, tmethod)
    testn += 1


# Create individual tests.  You may restrict which tests are generated
# by replacing the sequences in the ``for`` statements.  For instance:
testn = 0
for type_ in type_info:  # for type_ in ['string']:
    for op in operators:  # for op in ['!=']:
        for extracond in extra_conditions:  # for extracond in ['']:
            add_test_method(type_, op, extracond)

for type_ in ["float32", "float64"]:
    for func in func_info:  # i for func in ['log10']:
        for op in operators:
            add_test_method(type_, op, func=func)

# Base classes for non-indexed queries.
NX_BLOCK_SIZE1 = 128  # from ``interpreter.c`` in Numexpr
NX_BLOCK_SIZE2 = 8  # from ``interpreter.c`` in Numexpr


class SmallNITableMixin:
    nrows = row_period * 2
    assert NX_BLOCK_SIZE2 < nrows < NX_BLOCK_SIZE1
    assert nrows % NX_BLOCK_SIZE2 != 0  # to have some residual rows


class BigNITableMixin:
    nrows = row_period * 3
    assert nrows > NX_BLOCK_SIZE1 + NX_BLOCK_SIZE2
    assert nrows % NX_BLOCK_SIZE1 != 0
    assert nrows % NX_BLOCK_SIZE2 != 0  # to have some residual rows


# Parameters for non-indexed queries.
table_sizes = ["Small", "Big"]
heavy_table_sizes = frozenset(["Big"])
table_ndims = ["Scalar"]  # to enable multidimensional testing, include 'MD'

# Non-indexed queries: ``[SB][SM]TDTestCase``, where:
#
# 1. S is for small and B is for big size table.
#    Sizes are listed in `table_sizes`.
# 2. S is for scalar and M for multidimensional columns.
#    Dimensionalities are listed in `table_ndims`.


def niclassdata():
    for size in table_sizes:
        heavy = size in heavy_table_sizes
        for ndim in table_ndims:
            classname = f"{size[0]}{ndim[0]}TDTestCase"
            cbasenames = (
                f"{size}NITableMixin",
                f"{ndim}TableMixin",
                "TableDataTestCase",
            )
            classdict = dict(heavy=heavy)
            yield (classname, cbasenames, classdict)


# Base classes for the different type index.
class UltraLightITableMixin:
    kind = "ultralight"


class LightITableMixin:
    kind = "light"


class MediumITableMixin:
    kind = "medium"


class FullITableMixin:
    kind = "full"


# Base classes for indexed queries.


class SmallSTableMixin:
    nrows = 50


class MediumSTableMixin:
    nrows = 100


class BigSTableMixin:
    nrows = 500


# Parameters for indexed queries.
ckinds = ["UltraLight", "Light", "Medium", "Full"]
itable_sizes = ["Small", "Medium", "Big"]
heavy_itable_sizes = frozenset(["Medium", "Big"])
itable_optvalues = [0, 1, 3, 7, 9]
heavy_itable_optvalues = frozenset([0, 1, 7, 9])

# Indexed queries: ``[SMB]I[ulmf]O[01379]TDTestCase``, where:
#
# 1. S is for small, M for medium and B for big size table.
#    Sizes are listed in `itable_sizes`.
# 2. U is for 'ultraLight', L for 'light', M for 'medium', F for 'Full' indexes
#    Index types are listed in `ckinds`.
# 3. 0 to 9 is the desired index optimization level.
#    Optimizations are listed in `itable_optvalues`.


def iclassdata():
    for ckind in ckinds:
        for size in itable_sizes:
            for optlevel in itable_optvalues:
                heavy = (
                    optlevel in heavy_itable_optvalues
                    or size in heavy_itable_sizes
                )
                classname = "%sI%sO%dTDTestCase" % (
                    size[0],
                    ckind[0],
                    optlevel,
                )
                cbasenames = (
                    "%sSTableMixin" % size,
                    "%sITableMixin" % ckind,
                    "ScalarTableMixin",
                    "TableDataTestCase",
                )
                classdict = dict(heavy=heavy, optlevel=optlevel, indexed=True)
                yield (classname, cbasenames, classdict)


# Create test classes.
for cdatafunc in [niclassdata, iclassdata]:
    for cname, cbasenames, cdict in cdatafunc():
        cbases = tuple(eval(cbase) for cbase in cbasenames)
        class_ = type(cname, cbases, cdict)
        exec("%s = class_" % cname)


# Test cases on query usage
# -------------------------
class BaseTableUsageTestCase(BaseTableQueryTestCase):
    nrows = row_period


_gvar = None
"""Use this when a global variable is needed."""


class ScalarTableUsageTestCase(ScalarTableMixin, BaseTableUsageTestCase):
    """Test case for query usage on scalar tables.

    This also tests for most usage errors and situations.

    """

    def test_empty_condition(self):
        """Using an empty condition."""

        self.assertRaises(SyntaxError, self.table.where, "")

    def test_syntax_error(self):
        """Using a condition with a syntax error."""

        self.assertRaises(SyntaxError, self.table.where, "foo bar")

    def test_unsupported_object(self):
        """Using a condition with an unsupported object."""

        self.assertRaises((TypeError, ValueError), self.table.where, "[]")
        self.assertRaises(TypeError, self.table.where, "obj", {"obj": {}})
        self.assertRaises(
            (TypeError, ValueError), self.table.where, "c_bool < []"
        )

    def test_unsupported_syntax(self):
        """Using a condition with unsupported syntax."""

        self.assertRaises(
            (TypeError, ValueError), self.table.where, "c_bool[0]"
        )
        self.assertRaises(TypeError, self.table.where, "c_bool()")
        self.assertRaises(NameError, self.table.where, "c_bool.__init__")

    def test_no_column(self):
        """Using a condition with no participating columns."""

        self.assertRaises(ValueError, self.table.where, "True")

    def test_foreign_column(self):
        """Using a condition with a column from other table."""

        table2 = self.h5file.create_table("/", "other", self.tableDescription)
        self.assertRaises(
            ValueError,
            self.table.where,
            "c_int32_a + c_int32_b > 0",
            {
                "c_int32_a": self.table.cols.c_int32,
                "c_int32_b": table2.cols.c_int32,
            },
        )

    def test_unsupported_op(self):
        """Using a condition with unsupported operations on types."""

        NIE = NotImplementedError
        self.assertRaises(NIE, self.table.where, "c_complex128 > 0j")
        self.assertRaises(NIE, self.table.where, 'c_string + b"a" > b"abc"')

    def test_not_boolean(self):
        """Using a non-boolean condition."""

        self.assertRaises(TypeError, self.table.where, "c_int32")

    def test_nested_col(self):
        """Using a condition with nested columns."""

        self.assertRaises(TypeError, self.table.where, "c_nested")

    def test_implicit_col(self):
        """Using implicit column names in conditions."""

        # If implicit columns didn't work, a ``NameError`` would be raised.
        self.assertRaises(TypeError, self.table.where, "c_int32")
        # If overriding didn't work, no exception would be raised.
        self.assertRaises(
            TypeError,
            self.table.where,
            "c_bool",
            {"c_bool": self.table.cols.c_int32},
        )
        # External variables do not override implicit columns.

        def where_with_locals():
            c_int32 = self.table.cols.c_bool  # this wouldn't cause an error
            self.assertIsNotNone(c_int32)
            self.table.where("c_int32")

        self.assertRaises(TypeError, where_with_locals)

    def test_condition_vars(self):
        """Using condition variables in conditions."""

        # If condition variables didn't work, a ``NameError`` would be raised.
        self.assertRaises(
            NotImplementedError,
            self.table.where,
            "c_string > bound",
            {"bound": 0},
        )

        def where_with_locals():
            bound = "foo"  # this wouldn't cause an error
            # silence pyflakes warnings
            self.assertIsInstance(bound, str)
            self.table.where("c_string > bound", {"bound": 0})

        self.assertRaises(NotImplementedError, where_with_locals)

        def where_with_globals():
            global _gvar
            _gvar = "foo"  # this wouldn't cause an error
            # silence pyflakes warnings
            self.assertIsInstance(_gvar, str)
            try:
                self.table.where("c_string > _gvar", {"_gvar": 0})
            finally:
                del _gvar  # to keep global namespace clean

        self.assertRaises(NotImplementedError, where_with_globals)

    def test_scopes(self):
        """Looking up different scopes for variables."""

        # Make sure the variable is not implicit.
        self.assertRaises(NameError, self.table.where, "col")

        # First scope: dictionary of condition variables.
        self.assertRaises(
            TypeError,
            self.table.where,
            "col",
            {"col": self.table.cols.c_int32},
        )

        # Second scope: local variables.
        def where_whith_locals():
            col = self.table.cols.c_int32
            self.assertIsNotNone(col)
            self.table.where("col")

        self.assertRaises(TypeError, where_whith_locals)

        # Third scope: global variables.
        def where_with_globals():
            global _gvar
            _gvar = self.table.cols.c_int32
            # silence pyflakes warnings
            self.assertIsNotNone(_gvar)
            try:
                self.table.where("_gvar")
            finally:
                del _gvar  # to keep global namespace clean

        self.assertRaises(TypeError, where_with_globals)


class MDTableUsageTestCase(MDTableMixin, BaseTableUsageTestCase):
    """Test case for query usage on multidimensional tables."""

    def test(self):
        """Using a condition on a multidimensional table."""

        # Easy: queries on multidimensional tables are not implemented yet!
        self.assertRaises(NotImplementedError, self.table.where, "c_bool")


class IndexedTableUsage(ScalarTableMixin, BaseTableUsageTestCase):
    """Test case for query usage on indexed tables.

    Indexing could be used in more cases, but it is expected to kick in
    at least in the cases tested here.

    """

    nrows = 50
    indexed = True

    def setUp(self):
        super().setUp()
        self.table.cols.c_bool.create_index(_blocksizes=small_blocksizes)
        self.table.cols.c_int32.create_index(_blocksizes=small_blocksizes)
        self.will_query_use_indexing = self.table.will_query_use_indexing
        self.compileCondition = self.table._compile_condition
        self.requiredExprVars = self.table._required_expr_vars
        usable_idxs = set()
        for expr in self.idx_expr:
            idxvar = expr[0]
            if idxvar not in usable_idxs:
                usable_idxs.add(idxvar)
        self.usable_idxs = frozenset(usable_idxs)

    def test(self):
        for condition in self.conditions:
            c_usable_idxs = self.will_query_use_indexing(condition, {})
            self.assertEqual(
                c_usable_idxs,
                self.usable_idxs,
                f"\nQuery with condition: ``{condition}``\n"
                f"Computed usable indexes are: "
                f"``{c_usable_idxs}``\nand should be: "
                f"``{self.usable_idxs}``",
            )
            condvars = self.requiredExprVars(condition, None)
            compiled = self.compileCondition(condition, condvars)
            c_idx_expr = compiled.index_expressions
            self.assertEqual(
                c_idx_expr,
                self.idx_expr,
                f"\nWrong index expression in condition:\n"
                f"``{condition}``\nCompiled index expression is:"
                f"\n``{c_idx_expr}``\nand should be:\n"
                f"``{self.idx_expr}``",
            )
            c_str_expr = compiled.string_expression
            self.assertEqual(
                c_str_expr,
                self.str_expr,
                f"\nWrong index operations in condition:\n"
                f"``{condition}``\nComputed index operations are:"
                f"\n``{c_str_expr}``\nand should be:\n"
                f"``{self.str_expr}``",
            )
            common.verbosePrint(
                f"* Query with condition ``{condition}`` will use variables "
                f"``{compiled.index_variables}`` for indexing."
            )


class IndexedTableUsage1(IndexedTableUsage):
    conditions = [
        "(c_int32 > 0)",
        "(c_int32 > 0) & (c_extra > 0)",
        "(c_int32 > 0) & ((~c_bool) | (c_extra > 0))",
        "(c_int32 > 0) & ((c_extra < 3) & (c_extra > 0))",
    ]
    idx_expr = [("c_int32", ("gt",), (0,))]
    str_expr = "e0"


class IndexedTableUsage2(IndexedTableUsage):
    conditions = [
        "(c_int32 > 0) & (c_int32 < 5)",
        "(c_int32 > 0) & (c_int32 < 5) & (c_extra > 0)",
        "(c_int32 > 0) & (c_int32 < 5) & ((c_bool == True) | (c_extra > 0))",
        "(c_int32 > 0) & (c_int32 < 5) & ((c_extra > 0) | (c_bool == True))",
    ]
    idx_expr = [("c_int32", ("gt", "lt"), (0, 5))]
    str_expr = "e0"


class IndexedTableUsage3(IndexedTableUsage):
    conditions = [
        "(c_bool == True)",
        "(c_bool == True) & (c_extra > 0)",
        "(c_extra > 0) & (c_bool == True)",
        "((c_extra > 0) & (c_extra < 4)) & (c_bool == True)",
        "(c_bool == True) & ((c_extra > 0) & (c_extra < 4))",
    ]
    idx_expr = [("c_bool", ("eq",), (True,))]
    str_expr = "e0"


class IndexedTableUsage4(IndexedTableUsage):
    conditions = [
        "((c_int32 > 0) & (c_bool == True)) & (c_extra > 0)",
        "((c_int32 > 0) & (c_bool == True)) & ((c_extra > 0)"
        + " & (c_extra < 4))",
    ]
    idx_expr = [
        ("c_int32", ("gt",), (0,)),
        ("c_bool", ("eq",), (True,)),
    ]
    str_expr = "(e0 & e1)"


class IndexedTableUsage5(IndexedTableUsage):
    conditions = [
        "(c_int32 >= 1) & (c_int32 < 2) & (c_bool == True)",
        "(c_int32 >= 1) & (c_int32 < 2) & (c_bool == True)"
        + " & (c_extra > 0)",
    ]
    idx_expr = [
        ("c_int32", ("ge", "lt"), (1, 2)),
        ("c_bool", ("eq",), (True,)),
    ]
    str_expr = "(e0 & e1)"


class IndexedTableUsage6(IndexedTableUsage):
    conditions = [
        "(c_int32 >= 1) & (c_int32 < 2) & (c_int32 > 0) & (c_int32 < 5)",
        "(c_int32 >= 1) & (c_int32 < 2) & (c_int32 > 0) & (c_int32 < 5)"
        + " & (c_extra > 0)",
    ]
    idx_expr = [
        ("c_int32", ("ge", "lt"), (1, 2)),
        ("c_int32", ("gt",), (0,)),
        ("c_int32", ("lt",), (5,)),
    ]
    str_expr = "((e0 & e1) & e2)"


class IndexedTableUsage7(IndexedTableUsage):
    conditions = [
        "(c_int32 >= 1) & (c_int32 < 2) & ((c_int32 > 0) & (c_int32 < 5))",
        "((c_int32 >= 1) & (c_int32 < 2)) & ((c_int32 > 0) & (c_int32 < 5))",
        "((c_int32 >= 1) & (c_int32 < 2)) & ((c_int32 > 0) & (c_int32 < 5))"
        + " & (c_extra > 0)",
    ]
    idx_expr = [
        ("c_int32", ("ge", "lt"), (1, 2)),
        ("c_int32", ("gt", "lt"), (0, 5)),
    ]
    str_expr = "(e0 & e1)"


class IndexedTableUsage8(IndexedTableUsage):
    conditions = [
        "(c_extra > 0) & ((c_int32 > 0) & (c_int32 < 5))",
    ]
    idx_expr = [
        ("c_int32", ("gt", "lt"), (0, 5)),
    ]
    str_expr = "e0"


class IndexedTableUsage9(IndexedTableUsage):
    conditions = [
        "(c_extra > 0) & (c_int32 > 0) & (c_int32 < 5)",
        "((c_extra > 0) & (c_int32 > 0)) & (c_int32 < 5)",
        "(c_extra > 0) & (c_int32 > 0) & (c_int32 < 5) & (c_extra > 3)",
    ]
    idx_expr = [("c_int32", ("gt",), (0,)), ("c_int32", ("lt",), (5,))]
    str_expr = "(e0 & e1)"


class IndexedTableUsage10(IndexedTableUsage):
    conditions = [
        "(c_int32 < 5) & (c_extra > 0) & (c_bool == True)",
        "(c_int32 < 5) & (c_extra > 2) & c_bool",
        "(c_int32 < 5) & (c_bool == True) & (c_extra > 0) & (c_extra < 4)",
        "(c_int32 < 5) & (c_extra > 0) & (c_bool == True) & (c_extra < 4)",
    ]
    idx_expr = [("c_int32", ("lt",), (5,)), ("c_bool", ("eq",), (True,))]
    str_expr = "(e0 & e1)"


class IndexedTableUsage11(IndexedTableUsage):
    """Complex operations are not eligible for indexing."""

    conditions = [
        "sin(c_int32) > 0",
        "(c_int32 * 2.4) > 0",
        "(c_int32 + c_int32) > 0",
        "c_int32**2 > 0",
    ]
    idx_expr = []
    str_expr = ""


class IndexedTableUsage12(IndexedTableUsage):
    conditions = [
        "~c_bool",
        "~(c_bool)",
        "~c_bool & (c_extra > 0)",
        "~(c_bool) & (c_extra > 0)",
    ]
    idx_expr = [("c_bool", ("eq",), (False,))]
    str_expr = "e0"


class IndexedTableUsage13(IndexedTableUsage):
    conditions = [
        "~(c_bool == True)",
        "~((c_bool == True))",
        "~(c_bool == True) & (c_extra > 0)",
        "~((c_bool == True)) & (c_extra > 0)",
    ]
    idx_expr = [("c_bool", ("eq",), (False,))]
    str_expr = "e0"


class IndexedTableUsage14(IndexedTableUsage):
    conditions = [
        "~(c_int32 > 0)",
        "~((c_int32 > 0)) & (c_extra > 0)",
        "~(c_int32 > 0) & ((~c_bool) | (c_extra > 0))",
        "~(c_int32 > 0) & ((c_extra < 3) & (c_extra > 0))",
    ]
    idx_expr = [("c_int32", ("le",), (0,))]
    str_expr = "e0"


class IndexedTableUsage15(IndexedTableUsage):
    conditions = [
        "(~(c_int32 > 0) | ~c_bool)",
        "(~(c_int32 > 0) | ~(c_bool)) & (c_extra > 0)",
        "(~(c_int32 > 0) | ~(c_bool == True)) & ((c_extra > 0)"
        + " & (c_extra < 4))",
    ]
    idx_expr = [
        ("c_int32", ("le",), (0,)),
        ("c_bool", ("eq",), (False,)),
    ]
    str_expr = "(e0 | e1)"


class IndexedTableUsage16(IndexedTableUsage):
    conditions = [
        "(~(c_int32 > 0) & ~(c_int32 < 2))",
        "(~(c_int32 > 0) & ~(c_int32 < 2)) & (c_extra > 0)",
        "(~(c_int32 > 0) & ~(c_int32 < 2)) & ((c_extra > 0)"
        + " & (c_extra < 4))",
    ]
    idx_expr = [
        ("c_int32", ("le",), (0,)),
        ("c_int32", ("ge",), (2,)),
    ]
    str_expr = "(e0 & e1)"


class IndexedTableUsage17(IndexedTableUsage):
    conditions = [
        "(~(c_int32 > 0) & ~(c_int32 < 2))",
        "(~(c_int32 > 0) & ~(c_int32 < 2)) & (c_extra > 0)",
        "(~(c_int32 > 0) & ~(c_int32 < 2)) & ((c_extra > 0)"
        + " & (c_extra < 4))",
    ]
    idx_expr = [
        ("c_int32", ("le",), (0,)),
        ("c_int32", ("ge",), (2,)),
    ]
    str_expr = "(e0 & e1)"


# Negations of complex conditions are not supported yet


class IndexedTableUsage18(IndexedTableUsage):
    conditions = [
        "~((c_int32 > 0) & (c_bool))",
        "~((c_int32 > 0) & (c_bool)) & (c_extra > 0)",
        "~((c_int32 > 0) & (c_bool)) & ((c_extra > 0)" + " & (c_extra < 4))",
    ]
    idx_expr = []
    str_expr = ""


class IndexedTableUsage19(IndexedTableUsage):
    conditions = [
        "~((c_int32 > 0) & (c_bool)) & ((c_bool == False)"
        + " & (c_extra < 4))",
    ]
    idx_expr = [
        ("c_bool", ("eq",), (False,)),
    ]
    str_expr = "e0"


class IndexedTableUsage20(IndexedTableUsage):
    conditions = [
        "((c_int32 > 0) & ~(c_bool))",
        "((c_int32 > 0) & ~(c_bool)) & (c_extra > 0)",
        "((c_int32 > 0) & ~(c_bool == True)) & ((c_extra > 0) & (c_extra < 4))",
    ]
    idx_expr = [
        ("c_int32", ("gt",), (0,)),
        ("c_bool", ("eq",), (False,)),
    ]
    str_expr = "(e0 & e1)"


class IndexedTableUsage21(IndexedTableUsage):
    conditions = [
        "(~(c_int32 > 0) & (c_bool))",
        "(~(c_int32 > 0) & (c_bool)) & (c_extra > 0)",
        "(~(c_int32 > 0) & (c_bool == True)) & ((c_extra > 0)"
        + " & (c_extra < 4))",
    ]
    idx_expr = [
        ("c_int32", ("le",), (0,)),
        ("c_bool", ("eq",), (True,)),
    ]
    str_expr = "(e0 & e1)"


class IndexedTableUsage22(IndexedTableUsage):
    conditions = [
        "~((c_int32 >= 1) & (c_int32 < 2)) & ~(c_bool == True)",
        "~(c_bool == True) & (c_extra > 0)",
        "~((c_int32 >= 1) & (c_int32 < 2)) & (~(c_bool == True)"
        + " & (c_extra > 0))",
    ]
    idx_expr = [
        ("c_bool", ("eq",), (False,)),
    ]
    str_expr = "e0"


class IndexedTableUsage23(IndexedTableUsage):
    conditions = [
        "c_int32 != 1",
        "c_bool != False",
        "~(c_int32 != 1)",
        "~(c_bool != False)",
        "(c_int32 != 1) & (c_extra != 2)",
    ]
    idx_expr = []
    str_expr = ""


class IndexedTableUsage24(IndexedTableUsage):
    conditions = [
        "c_bool",
        "c_bool == True",
        "True == c_bool",
        "~(~c_bool)",
        "~~c_bool",
        "~~~~c_bool",
        "~(~c_bool) & (c_extra != 2)",
    ]
    idx_expr = [
        ("c_bool", ("eq",), (True,)),
    ]
    str_expr = "e0"


class IndexedTableUsage25(IndexedTableUsage):
    conditions = [
        "~c_bool",
        "c_bool == False",
        "False == c_bool",
        "~(c_bool)",
        "~((c_bool))",
        "~~~c_bool",
        "~~(~c_bool) & (c_extra != 2)",
    ]
    idx_expr = [
        ("c_bool", ("eq",), (False,)),
    ]
    str_expr = "e0"


class IndexedTableUsage26(IndexedTableUsage):
    conditions = [
        "c_bool != True",
        "True != c_bool",
        "c_bool != False",
        "False != c_bool",
    ]
    idx_expr = []
    str_expr = ""


class IndexedTableUsage27(IndexedTableUsage):
    conditions = [
        "(c_int32 == 3) | c_bool | (c_int32 == 5)",
        "(((c_int32 == 3) | (c_bool == True)) | (c_int32 == 5))"
        + " & (c_extra > 0)",
    ]
    idx_expr = [
        ("c_int32", ("eq",), (3,)),
        ("c_bool", ("eq",), (True,)),
        ("c_int32", ("eq",), (5,)),
    ]
    str_expr = "((e0 | e1) | e2)"


class IndexedTableUsage28(IndexedTableUsage):
    conditions = [
        "((c_int32 == 3) | c_bool) & (c_int32 == 5)",
        "(((c_int32 == 3) | (c_bool == True)) & (c_int32 == 5))"
        + " & (c_extra > 0)",
    ]
    idx_expr = [
        ("c_int32", ("eq",), (3,)),
        ("c_bool", ("eq",), (True,)),
        ("c_int32", ("eq",), (5,)),
    ]
    str_expr = "((e0 | e1) & e2)"


class IndexedTableUsage29(IndexedTableUsage):
    conditions = [
        "(c_int32 == 3) | ((c_int32 == 4) & (c_int32 == 5))",
        "((c_int32 == 3) | ((c_int32 == 4) & (c_int32 == 5)))"
        + " & (c_extra > 0)",
    ]
    idx_expr = [
        ("c_int32", ("eq",), (4,)),
        ("c_int32", ("eq",), (5,)),
        ("c_int32", ("eq",), (3,)),
    ]
    str_expr = "((e0 & e1) | e2)"


class IndexedTableUsage30(IndexedTableUsage):
    conditions = [
        "((c_int32 == 3) | (c_int32 == 4)) & (c_int32 == 5)",
        "((c_int32 == 3) | (c_int32 == 4)) & (c_int32 == 5)"
        + " & (c_extra > 0)",
    ]
    idx_expr = [
        ("c_int32", ("eq",), (3,)),
        ("c_int32", ("eq",), (4,)),
        ("c_int32", ("eq",), (5,)),
    ]
    str_expr = "((e0 | e1) & e2)"


class IndexedTableUsage31(IndexedTableUsage):
    conditions = [
        "(c_extra > 0) & ((c_extra < 4) & (c_bool == True))",
        "(c_extra > 0) & ((c_bool == True) & (c_extra < 5))",
        "((c_int32 > 0) | (c_extra > 0)) & (c_bool == True)",
    ]
    idx_expr = [
        ("c_bool", ("eq",), (True,)),
    ]
    str_expr = "e0"


class IndexedTableUsage32(IndexedTableUsage):
    conditions = [
        "(c_int32 < 5) & (c_extra > 0) & (c_bool == True) | (c_extra < 4)",
    ]
    idx_expr = []
    str_expr = ""


# Main part
# ---------
def suite():
    """Return a test suite consisting of all the test cases in the module."""

    testSuite = common.unittest.TestSuite()

    cdatafuncs = [niclassdata]  # non-indexing data tests
    cdatafuncs.append(iclassdata)  # indexing data tests

    heavy = common.heavy
    # Choose which tests to run in classes with autogenerated tests.
    if heavy:
        autoprefix = "test"  # all tests
    else:
        autoprefix = "test_l"  # only light tests

    niter = 1
    for i in range(niter):
        # Tests on query data.
        for cdatafunc in cdatafuncs:
            for cdata in cdatafunc():
                class_ = eval(cdata[0])
                if heavy or not class_.heavy:
                    suite_ = common.make_suite(class_, prefix=autoprefix)
                    testSuite.addTest(suite_)
        # Tests on query usage.
        testSuite.addTest(common.make_suite(ScalarTableUsageTestCase))
        testSuite.addTest(common.make_suite(MDTableUsageTestCase))
        testSuite.addTest(common.make_suite(IndexedTableUsage1))
        testSuite.addTest(common.make_suite(IndexedTableUsage2))
        testSuite.addTest(common.make_suite(IndexedTableUsage3))
        testSuite.addTest(common.make_suite(IndexedTableUsage4))
        testSuite.addTest(common.make_suite(IndexedTableUsage5))
        testSuite.addTest(common.make_suite(IndexedTableUsage6))
        testSuite.addTest(common.make_suite(IndexedTableUsage7))
        testSuite.addTest(common.make_suite(IndexedTableUsage8))
        testSuite.addTest(common.make_suite(IndexedTableUsage9))
        testSuite.addTest(common.make_suite(IndexedTableUsage10))
        testSuite.addTest(common.make_suite(IndexedTableUsage11))
        testSuite.addTest(common.make_suite(IndexedTableUsage12))
        testSuite.addTest(common.make_suite(IndexedTableUsage13))
        testSuite.addTest(common.make_suite(IndexedTableUsage14))
        testSuite.addTest(common.make_suite(IndexedTableUsage15))
        testSuite.addTest(common.make_suite(IndexedTableUsage16))
        testSuite.addTest(common.make_suite(IndexedTableUsage17))
        testSuite.addTest(common.make_suite(IndexedTableUsage18))
        testSuite.addTest(common.make_suite(IndexedTableUsage19))
        testSuite.addTest(common.make_suite(IndexedTableUsage20))
        testSuite.addTest(common.make_suite(IndexedTableUsage21))
        testSuite.addTest(common.make_suite(IndexedTableUsage22))
        testSuite.addTest(common.make_suite(IndexedTableUsage23))
        testSuite.addTest(common.make_suite(IndexedTableUsage24))
        testSuite.addTest(common.make_suite(IndexedTableUsage25))
        testSuite.addTest(common.make_suite(IndexedTableUsage26))
        testSuite.addTest(common.make_suite(IndexedTableUsage27))
        testSuite.addTest(common.make_suite(IndexedTableUsage28))
        testSuite.addTest(common.make_suite(IndexedTableUsage29))
        testSuite.addTest(common.make_suite(IndexedTableUsage30))
        testSuite.addTest(common.make_suite(IndexedTableUsage31))
        testSuite.addTest(common.make_suite(IndexedTableUsage32))

    return testSuite


if __name__ == "__main__":
    common.parse_argv(sys.argv)
    common.print_versions()
    common.unittest.main(defaultTest="suite")
