"""Test module for evaluating expressions under PyTables."""

import numpy as np
from numpy import testing as npt

import tables as tb
from tables.tests import common

# An example of record


class Record(tb.IsDescription):
    colInt32 = tb.Int32Col()
    colInt64 = tb.Int64Col()
    colFloat32 = tb.Float32Col()
    colFloat64 = tb.Float64Col()
    colComplex = tb.ComplexCol(itemsize=16)


# Helper functions
def get_sliced_vars(npvars, start, stop, step):
    npvars_ = {}
    for name, var in npvars.items():
        if hasattr(var, "__len__"):
            npvars_[name] = var[start:stop:step]
        else:
            npvars_[name] = var
    return npvars_


def get_sliced_vars2(npvars, start, stop, step, shape, maindim):
    npvars_ = {}
    slices = [slice(None) for dim in shape]
    slices[maindim] = slice(start, stop, step)
    for name, var in npvars.items():
        npvars_[name] = var.__getitem__(tuple(slices))
    return npvars_


# Basic tests
class ExprTestCase(common.TempFileMixin, common.PyTablesTestCase):

    # The shape for the variables in expressions
    shape = (10, 20)

    def setUp(self):
        super().setUp()

        # The expression
        self.expr = "2 * a*b + c"
        # Define the NumPy variables to be used in expression
        N = np.prod(self.shape)
        self.a = a = np.arange(0, N, dtype="int32").reshape(self.shape)
        self.b = b = np.arange(N, 2 * N, dtype="int64").reshape(self.shape)
        self.c = c = np.arange(2 * N, 3 * N, dtype="int32").reshape(self.shape)
        self.r1 = r1 = np.empty(N, dtype="int64").reshape(self.shape)
        self.npvars = {
            "a": a,
            "b": b,
            "c": c,
        }
        # Define other variables, if needed
        root = self.h5file.root
        if self.kind == "Array":
            self.a = self.h5file.create_array(root, "a", a)
            self.b = self.h5file.create_array(root, "b", b)
            self.c = self.h5file.create_array(root, "c", c)
            self.r1 = self.h5file.create_array(root, "r1", r1)
        elif self.kind == "CArray":
            self.a = self.h5file.create_carray(
                root, "a", atom=tb.Atom.from_dtype(a.dtype), shape=self.shape
            )
            self.b = self.h5file.create_carray(
                root, "b", atom=tb.Atom.from_dtype(b.dtype), shape=self.shape
            )
            self.c = self.h5file.create_carray(
                root, "c", atom=tb.Atom.from_dtype(c.dtype), shape=self.shape
            )
            self.r1 = self.h5file.create_carray(
                root, "r1", atom=tb.Atom.from_dtype(r1.dtype), shape=self.shape
            )
            self.a[:] = a
            self.b[:] = b
            self.c[:] = c
        elif self.kind == "EArray":
            shape = list(self.shape)
            shape[0] = 0
            self.a = self.h5file.create_earray(
                root, "a", atom=tb.Atom.from_dtype(a.dtype), shape=shape
            )
            self.b = self.h5file.create_earray(
                root, "b", atom=tb.Atom.from_dtype(b.dtype), shape=shape
            )
            self.c = self.h5file.create_earray(
                root, "c", atom=tb.Atom.from_dtype(c.dtype), shape=shape
            )
            self.r1 = self.h5file.create_earray(
                root, "r1", atom=tb.Atom.from_dtype(r1.dtype), shape=shape
            )
            self.a.append(a)
            self.b.append(b)
            self.c.append(c)
            self.r1.append(r1)  # Fill with uninitialized values
        elif self.kind == "Column":
            ra = np.rec.fromarrays(
                [a, b, c, r1],
                dtype="%si4,%si8,%si4,%si8" % ((self.shape[1:],) * 4),
            )
            t = self.h5file.create_table(root, "t", ra)
            self.a = t.cols.f0
            self.b = t.cols.f1
            self.c = t.cols.f2
            self.d = t.cols.f3
        self.vars = {
            "a": self.a,
            "b": self.b,
            "c": self.c,
        }

    def test00_simple(self):
        """Checking that expression is correctly evaluated."""

        expr = tb.Expr(self.expr, self.vars)
        r1 = expr.eval()
        r2 = eval(self.expr, self.npvars)
        if common.verbose:
            print("Computed expression:", repr(r1))
            print("Should look like:", repr(r2))
        self.assertTrue(
            common.areArraysEqual(r1, r2),
            "Evaluate is returning a wrong value.",
        )

    def test01_out(self):
        """Checking that expression is correctly evaluated (`out` param)"""

        expr = tb.Expr(self.expr, self.vars)
        expr.set_output(self.r1)
        r1 = expr.eval()
        if self.kind != "NumPy":
            r1 = r1[:]
        r2 = eval(self.expr, self.npvars)
        if common.verbose:
            print("Computed expression:", repr(r1))
            print("Should look like:", repr(r2))
        self.assertTrue(
            common.areArraysEqual(r1, r2),
            "Evaluate is returning a wrong value.",
        )

    def test02_out(self):
        """Checking that expression is correctly evaluated when slice is
        outside of data samples (`out` param)"""
        expr = tb.Expr(self.expr, self.vars)
        # maybe it's better to use the leading dimension instead?
        maxshape = max(self.shape)
        start, stop, step = (maxshape + 1, maxshape + 2, None)
        expr.set_inputs_range(start, stop, step)
        r1 = expr.eval()
        # create an empty array with the same dtype and shape
        zeros = np.zeros(shape=self.shape, dtype=r1.dtype)
        r2 = zeros[start:stop:step]
        self.assertListEqual(r1.tolist(), r2.tolist())
        if common.verbose:
            print("Computed expression:", repr(r1))
            print("Should look like:", repr(r2))
        self.assertTrue(
            common.areArraysEqual(r1, r2),
            "Evaluate is returning a wrong value.",
        )


class ExprNumPy(ExprTestCase):
    kind = "NumPy"


class ExprArray(ExprTestCase):
    kind = "Array"


class ExprCArray(ExprTestCase):
    kind = "CArray"


class ExprEArray(ExprTestCase):
    kind = "EArray"


class ExprColumn(ExprTestCase):
    kind = "Column"


# Test for mixed containers
class MixedContainersTestCase(common.TempFileMixin, common.PyTablesTestCase):

    def setUp(self):
        super().setUp()

        # The expression
        self.expr = "2 * a*b + c**2+d**2+e-f+g"

        # Create a directory in file for outputs
        root = self.h5file.root
        outs = self.h5file.create_group(root, "outs")

        # Define the NumPy variables to be used in expression
        N = np.prod(self.shape)

        # Initial values for variables
        a = np.arange(0, N, dtype="int32").reshape(self.shape)
        b = np.arange(N, 2 * N, dtype="int64").reshape(self.shape)
        c = np.arange(2 * N, 3 * N, dtype="int32").reshape(self.shape)
        d = np.arange(3 * N, 4 * N, dtype="int32").reshape(self.shape)
        e = np.arange(4 * N, 5 * N, dtype="int32").reshape(self.shape)
        self.f = f = int(3)  # a regular python type
        self.g = g = np.int16(2)  # a NumPy scalar type

        # Original values
        self.npvars = {"a": a, "b": b, "c": c, "d": d, "e": e, "f": f, "g": g}
        rnda = b.copy()

        # ndarray input and output
        self.a = a
        self.rnda = rnda

        # Array input and output
        self.b = self.h5file.create_array(root, "b", b)
        self.rarr = self.b.copy(outs)

        # CArray input and output
        self.c = self.h5file.create_carray(
            root, "c", atom=tb.Atom.from_dtype(c.dtype), shape=self.shape
        )
        self.c[:] = c
        self.rcarr = self.c.copy(outs)

        # EArray input and output
        eshape = list(self.shape)
        eshape[0] = 0
        self.d = self.h5file.create_earray(
            root, "d", atom=tb.Atom.from_dtype(d.dtype), shape=eshape
        )
        self.d.append(d)
        self.rearr = self.d.copy(outs)

        # Column input and output
        rtype = {}
        colshape = self.shape[1:]
        for i, col in enumerate((a, b, c, d, e, rnda)):
            rtype["f%d" % i] = tb.Col.from_sctype(col.dtype.type, colshape)
        t = self.h5file.create_table(root, "t", rtype)
        nrows = self.shape[0]
        row = t.row
        for nrow in range(nrows):
            for i, col in enumerate((a, b, c, d, e, rnda)):
                row["f%d" % i] = col[nrow]
            row.append()
        t.flush()
        self.e = t.cols.f4
        self.rcol = t.cols.f5
        # Input vars
        self.vars = {
            "a": self.a,
            "b": self.b,
            "c": self.c,
            "d": self.d,
            "e": self.e,
            "f": self.f,
            "g": self.g,
        }

    def test00a_simple(self):
        """Checking expressions with mixed objects."""

        expr = tb.Expr(self.expr, self.vars)
        r1 = expr.eval()
        r2 = eval(self.expr, self.npvars)
        if common.verbose:
            print("Computed expression:", repr(r1), r1.dtype)
            print("Should look like:", repr(r2), r2.dtype)

        self.assertTrue(
            common.areArraysEqual(r1, r2),
            "Evaluate is returning a wrong value.",
        )

    def test00b_simple_scalars(self):
        """Checking that scalars in expression evaluate correctly."""

        expr_str = "2 * f + g"
        expr = tb.Expr(expr_str, self.vars)
        r1 = expr.eval()
        r2 = eval(expr_str, self.npvars)
        if common.verbose:
            print("Computed expression:", repr(r1), r1.dtype)
            print("Should look like:", repr(r2), r2.dtype)
        msg = f"Evaluate is returning a wrong value: {expr_str}\n{r1=}\n{r2=}"
        self.assertEqual(r1.shape, r2.shape, msg=msg)
        # In something like 2 * np.in16(3) + np.int16(2) the result is still a
        # np.int16 in NumPy 2.0, so we shouldn't actually check dtype but just
        # the kind
        self.assertEqual(r1.dtype.kind, r2.dtype.kind, msg=msg)
        self.assertEqual(r1, r2, msg=msg)

    def test01a_out(self):
        """Checking expressions with mixed objects (`out` param)"""

        expr = tb.Expr(self.expr, self.vars)
        for r1 in self.rnda, self.rarr, self.rcarr, self.rearr, self.rcol:
            if common.verbose:
                print("Checking output container:", type(r1))
            expr.set_output(r1)
            r1 = expr.eval()
            if not isinstance(r1, type(self.rnda)):
                r1 = r1[:]
            r2 = eval(self.expr, self.npvars)
            if common.verbose:
                print("Computed expression:", repr(r1), r1.dtype)
                print("Should look like:", repr(r2), r2.dtype)
            self.assertTrue(
                common.areArraysEqual(r1, r2),
                "Evaluate is returning a wrong value.",
            )

    def test01b_out_scalars(self):
        """Checking expressions with mixed objects (`out` param, scalars)"""

        if len(self.shape) > 1:
            # This test is only meant for undimensional outputs
            return
        expr_str = "2 * f + g"
        expr = tb.Expr(expr_str, self.vars)
        for r1 in self.rnda, self.rarr, self.rcarr, self.rearr, self.rcol:
            if common.verbose:
                print("Checking output container:", type(r1))
            expr.set_output(r1)
            r1 = expr.eval()
            r1 = r1[()]  # convert a 0-dim array into a scalar
            r2 = eval(expr_str, self.npvars)
            if common.verbose:
                print("Computed expression:", repr(r1), r1.dtype)
                print("Should look like:", repr(r2), r2.dtype)
            msg = (
                f"Evaluate is returning a wrong value: "
                f"{expr_str}\n{r1=}\n{r2=}"
            )
            # On NumPy 2 type promotion is different so don't check type
            # strictly here
            self.assertTrue(
                common.areArraysEqual(r1, r2, check_type=False), msg=msg
            )
            self.assertEqual(r1.dtype.kind, r2.dtype.kind)

    def test02a_sss(self):
        """Checking mixed objects and start, stop, step (I)"""

        start, stop, step = (self.start, self.stop, 1)
        expr = tb.Expr(self.expr, self.vars)
        expr.set_inputs_range(start, stop, step)
        r1 = expr.eval()
        npvars = get_sliced_vars(self.npvars, start, stop, step)
        r2 = eval(self.expr, npvars)
        if common.verbose:
            print("Computed expression:", repr(r1), r1.dtype)
            print("Should look like:", repr(r2), r2.dtype)
        self.assertTrue(
            common.areArraysEqual(r1, r2),
            "Evaluate is returning a wrong value.",
        )

    def test02b_sss(self):
        """Checking mixed objects and start, stop, step (II)"""

        start, stop, step = (0, self.shape[0], self.step)
        expr = tb.Expr(self.expr, self.vars)
        expr.set_inputs_range(start, stop, step)
        r1 = expr.eval()
        npvars = get_sliced_vars(self.npvars, start, stop, step)
        r2 = eval(self.expr, npvars)
        if common.verbose:
            print("Computed expression:", repr(r1), r1.dtype)
            print("Should look like:", repr(r2), r2.dtype)
        self.assertTrue(
            common.areArraysEqual(r1, r2),
            "Evaluate is returning a wrong value.",
        )

    def test02c_sss(self):
        """Checking mixed objects and start, stop, step (III)"""

        start, stop, step = (self.start, self.stop, self.step)
        expr = tb.Expr(self.expr, self.vars)
        expr.set_inputs_range(start, stop, step)
        r1 = expr.eval()
        npvars = get_sliced_vars(self.npvars, start, stop, step)
        r2 = eval(self.expr, npvars)
        if common.verbose:
            print("Computed expression:", repr(r1), r1.dtype)
            print("Should look like:", repr(r2), r2.dtype)
        self.assertTrue(
            common.areArraysEqual(r1, r2),
            "Evaluate is returning a wrong value.",
        )

    def test03_sss(self):
        """Checking start, stop, step as numpy.int64."""

        start, stop, step = (
            np.int64(i) for i in (self.start, self.stop, self.step)
        )
        expr = tb.Expr(self.expr, self.vars)
        expr.set_inputs_range(start, stop, step)
        r1 = expr.eval()
        npvars = get_sliced_vars(self.npvars, start, stop, step)
        r2 = eval(self.expr, npvars)
        if common.verbose:
            print("Computed expression:", repr(r1), r1.dtype)
            print("Should look like:", repr(r2), r2.dtype)
        self.assertTrue(
            common.areArraysEqual(r1, r2),
            "Evaluate is returning a wrong value.",
        )


class MixedContainers0(MixedContainersTestCase):
    shape = (1,)
    start, stop, step = (0, 1, 1)


class MixedContainers1(MixedContainersTestCase):
    shape = (10,)
    start, stop, step = (3, 6, 2)


class MixedContainers2(MixedContainersTestCase):
    shape = (10, 5)
    start, stop, step = (2, 9, 3)


class MixedContainers3(MixedContainersTestCase):
    shape = (10, 3, 2)
    start, stop, step = (2, -1, 1)


# Test for unaligned objects
class UnalignedObject(common.PyTablesTestCase):

    def test00_simple(self):
        """Checking expressions with unaligned objects."""

        # Build unaligned arrays
        a0 = np.empty(10, dtype="int8")
        a1 = np.arange(10, dtype="int32")
        a2 = a1.copy()
        a3 = a2.copy()
        ra = np.rec.fromarrays([a0, a1, a2, a3])
        # The inputs
        a = ra["f1"]
        b = ra["f2"]
        self.assertEqual(a.flags.aligned, False)
        self.assertEqual(b.flags.aligned, False)
        # The expression
        sexpr = "2 * a + b"
        expr = tb.Expr(sexpr)
        r1 = expr.eval()
        r2 = eval(sexpr)
        if common.verbose:
            print("Computed expression:", repr(r1), r1.dtype)
            print("Should look like:", repr(r2), r2.dtype)
        self.assertTrue(
            common.areArraysEqual(r1, r2),
            "Evaluate is returning a wrong value.",
        )

    def test01_md(self):
        """Checking expressions with unaligned objects (MD version)"""

        # Build unaligned arrays
        a0 = np.empty((10, 4), dtype="int8")
        a1 = np.arange(10 * 4, dtype="int32").reshape(10, 4)
        a2 = a1.copy()
        a3 = a2.copy()
        ra = np.rec.fromarrays([a0, a1, a2, a3])
        # The inputs
        a = ra["f1"]
        b = ra["f2"]
        self.assertEqual(a.flags.aligned, False)
        self.assertEqual(b.flags.aligned, False)
        # The expression
        sexpr = "2 * a + b"
        expr = tb.Expr(sexpr)
        r1 = expr.eval()
        r2 = eval(sexpr)
        if common.verbose:
            print("Computed expression:", repr(r1), r1.dtype)
            print("Should look like:", repr(r2), r2.dtype)
        self.assertTrue(
            common.areArraysEqual(r1, r2),
            "Evaluate is returning a wrong value.",
        )


# Test for non-contiguous objects
class NonContiguousObject(common.PyTablesTestCase):

    def test00_simple(self):
        """Checking expressions with non-contiguous objects"""

        # Build non-contiguous arrays as inputs
        a = np.arange(10, dtype="int32")
        b = a[::2]
        a = b * 2
        self.assertEqual(b.flags.contiguous, False)
        self.assertEqual(b.flags.aligned, True)
        # The expression
        sexpr = "2 * a + b"
        expr = tb.Expr(sexpr)
        r1 = expr.eval()
        r2 = eval(sexpr)
        if common.verbose:
            print("Computed expression:", repr(r1), r1.dtype)
            print("Should look like:", repr(r2), r2.dtype)
        self.assertTrue(
            common.areArraysEqual(r1, r2),
            "Evaluate is returning a wrong value.",
        )

    def test01a_md(self):
        """Checking expressions with non-contiguous objects (MD version, I)"""

        # Build non-contiguous arrays
        a = np.arange(10 * 4, dtype="int32").reshape(10, 4)
        b = a[::2]
        a = b * 2
        self.assertEqual(b.flags.contiguous, False)
        self.assertEqual(b.flags.aligned, True)
        # The expression
        sexpr = "2 * a + b"
        expr = tb.Expr(sexpr)
        r1 = expr.eval()
        r2 = eval(sexpr)
        if common.verbose:
            print("Computed expression:", repr(r1), r1.dtype)
            print("Should look like:", repr(r2), r2.dtype)
        self.assertTrue(
            common.areArraysEqual(r1, r2),
            "Evaluate is returning a wrong value.",
        )

    def test01b_md(self):
        """Checking expressions with non-contiguous objects (MD version, II)"""

        # Build non-contiguous arrays
        a = np.arange(10 * 4, dtype="int32").reshape(10, 4)
        b = a[:, ::2]
        a = b * 2
        self.assertEqual(b.flags.contiguous, False)
        self.assertEqual(b.flags.aligned, True)
        # The expression
        sexpr = "2 * a + b"
        expr = tb.Expr(sexpr)
        r1 = expr.eval()
        r2 = eval(sexpr)
        if common.verbose:
            print("Computed expression:", repr(r1), r1.dtype)
            print("Should look like:", repr(r2), r2.dtype)
        self.assertTrue(
            common.areArraysEqual(r1, r2),
            "Evaluate is returning a wrong value.",
        )


# Test for errors
class ExprError(common.TempFileMixin, common.PyTablesTestCase):

    # The shape for the variables in expressions
    shape = (10,)

    def setUp(self):
        super().setUp()

        # Define the NumPy variables to be used in expression
        N = np.prod(self.shape)
        self.a = np.arange(N, dtype="int32").reshape(self.shape)
        self.b = np.arange(N, dtype="int64").reshape(self.shape)
        self.c = np.arange(N, dtype="int32").reshape(self.shape)
        self.r1 = np.empty(N, dtype="int64").reshape(self.shape)

    def _test00_shape(self):
        """Checking that inconsistent shapes are detected."""

        self.b = self.b.reshape(self.shape + (1,))
        expr = "a * b + c"
        vars_ = {
            "a": self.a,
            "b": self.b,
            "c": self.c,
        }
        expr = tb.Expr(expr, vars_)
        self.assertRaises(ValueError, expr.eval)

    def test02_uint64(self):
        """Checking that uint64 arrays in expression are detected."""

        self.b = self.b.view("uint64")
        expr = "a * b + c"
        vars_ = {
            "a": self.a,
            "b": self.b,
            "c": self.c,
        }
        self.assertRaises(NotImplementedError, tb.Expr, expr, vars_)

    def test03_table(self):
        """Checking that tables in expression are detected."""

        class Rec(tb.IsDescription):
            col1 = tb.Int32Col()
            col2 = tb.Int64Col()

        t = self.h5file.create_table("/", "a", Rec)
        expr = "a * b + c"
        vars_ = {
            "a": t,
            "b": self.b,
            "c": self.c,
        }
        self.assertRaises(TypeError, tb.Expr, expr, vars_)

    def test04_nestedcols(self):
        """Checking that nested cols in expression are detected."""

        class Nested(tb.IsDescription):
            col1 = tb.Int32Col()

            class col2(tb.IsDescription):
                col3 = tb.Int64Col()

        t = self.h5file.create_table("/", "a", Nested)
        expr = "a * b + c"
        # The next non-nested column should work
        a = t.cols.col2.col3
        vars_ = {
            "a": a,
            "b": self.b,
            "c": self.c,
        }
        expr = tb.Expr(expr, vars_)
        r1 = expr.eval()
        self.assertIsNotNone(r1)
        # But a nested column should not
        a = t.cols.col2
        vars_ = {
            "a": a,
            "b": self.b,
            "c": self.c,
        }
        self.assertRaises(TypeError, tb.Expr, expr, vars_)

    def test05_vlarray(self):
        """Checking that VLArrays in expression are detected."""

        vla = self.h5file.create_vlarray("/", "a", tb.Int32Col())
        expr = "a * b + c"
        vars_ = {
            "a": vla,
            "b": self.b,
            "c": self.c,
        }
        self.assertRaises(TypeError, tb.Expr, expr, vars_)


# Test for broadcasting arrays
class BroadcastTestCase(common.TempFileMixin, common.PyTablesTestCase):

    def test00_simple(self):
        """Checking broadcast in expression."""

        shapes = (self.shape1, self.shape2, self.shape3)
        # Build arrays with different shapes as inputs
        a = np.arange(np.prod(shapes[0]), dtype="i4").reshape(shapes[0])
        b = np.arange(np.prod(shapes[1]), dtype="i4").reshape(shapes[1])
        c = np.arange(np.prod(shapes[2]), dtype="i4").reshape(shapes[2])
        root = self.h5file.root
        if a.shape[0] > 0:
            a1 = self.h5file.create_array(root, "a1", a)
        else:
            a1 = self.h5file.create_earray(
                root, "a1", atom=tb.Int32Col(), shape=a.shape
            )
        self.assertIsNotNone(a1)
        b1 = self.h5file.create_array(root, "b1", b)
        self.assertIsNotNone(b1)
        c1 = self.h5file.create_array(root, "c1", c)
        self.assertIsNotNone(c1)
        # The expression
        expr = tb.Expr("2 * a1 + b1-c1")
        r1 = expr.eval()
        r2 = eval("2 * a + b-c")
        if common.verbose:
            print("Tested shapes:", self.shape1, self.shape2, self.shape3)
            print("Computed expression:", repr(r1), r1.dtype)
            print("Should look like:", repr(r2), r2.dtype)
        self.assertTrue(
            common.areArraysEqual(r1, r2),
            "Evaluate is returning a wrong value.",
        )


class Broadcast0(BroadcastTestCase):
    shape1 = (0, 3, 4)
    shape2 = (3, 4)
    shape3 = (4,)


class Broadcast1(BroadcastTestCase):
    shape1 = (2, 3, 4)
    shape2 = (3, 4)
    shape3 = (4,)


class Broadcast2(BroadcastTestCase):
    shape1 = (
        3,
        4,
    )
    shape2 = (3, 4)
    shape3 = (4,)


class Broadcast3(BroadcastTestCase):
    shape1 = (4,)
    shape2 = (3, 4)
    shape3 = (4,)


class Broadcast4(BroadcastTestCase):
    shape1 = (1,)
    shape2 = (3, 4)
    shape3 = (4,)


class Broadcast5(BroadcastTestCase):
    shape1 = (1,)
    shape2 = (3, 1)
    shape3 = (4,)


# Test for different length inputs
class DiffLengthTestCase(common.TempFileMixin, common.PyTablesTestCase):

    def test00_simple(self):
        """Checking different length inputs in expression."""

        shapes = (list(self.shape1), list(self.shape2), list(self.shape3))
        # Build arrays with different shapes as inputs
        a = np.arange(np.prod(shapes[0]), dtype="i4").reshape(shapes[0])
        b = np.arange(np.prod(shapes[1]), dtype="i4").reshape(shapes[1])
        c = np.arange(np.prod(shapes[2]), dtype="i4").reshape(shapes[2])
        # The expression
        expr = tb.Expr("2 * a + b-c")
        r1 = expr.eval()
        # Compute the minimum length for shapes
        maxdim = max([len(shape) for shape in shapes])
        minlen = min(
            [
                shape[0]
                for i, shape in enumerate(shapes)
                if len(shape) == maxdim
            ]
        )
        for i, shape in enumerate(shapes):
            if len(shape) == maxdim:
                shape[0] = minlen
        # Build arrays with the new shapes as inputs
        a = np.arange(np.prod(shapes[0]), dtype="i4").reshape(shapes[0])
        self.assertIsNotNone(a)
        b = np.arange(np.prod(shapes[1]), dtype="i4").reshape(shapes[1])
        self.assertIsNotNone(b)
        c = np.arange(np.prod(shapes[2]), dtype="i4").reshape(shapes[2])
        self.assertIsNotNone(c)
        r2 = eval("2 * a + b-c")
        if common.verbose:
            print("Tested shapes:", self.shape1, self.shape2, self.shape3)
            print("Computed expression:", repr(r1), r1.dtype)
            print("Should look like:", repr(r2), r2.dtype)
        self.assertTrue(
            common.areArraysEqual(r1, r2),
            "Evaluate is returning a wrong value.",
        )


class DiffLength0(DiffLengthTestCase):
    shape1 = (0,)
    shape2 = (10,)
    shape3 = (20,)


class DiffLength1(DiffLengthTestCase):
    shape1 = (3,)
    shape2 = (10,)
    shape3 = (20,)


class DiffLength2(DiffLengthTestCase):
    shape1 = (3, 4)
    shape2 = (2, 3, 4)
    shape3 = (4, 3, 4)


class DiffLength3(DiffLengthTestCase):
    shape1 = (1, 3, 4)
    shape2 = (2, 3, 4)
    shape3 = (4, 3, 4)


class DiffLength4(DiffLengthTestCase):
    shape1 = (0, 3, 4)
    shape2 = (2, 3, 4)
    shape3 = (4, 3, 4)


# Test for different type inputs
class TypesTestCase(common.TempFileMixin, common.PyTablesTestCase):

    def test00_bool(self):
        """Checking booleans in expression."""

        # Build arrays with different shapes as inputs
        a = np.array([True, False, True])
        b = np.array([False, True, False])
        root = self.h5file.root
        a1 = self.h5file.create_array(root, "a1", a)
        self.assertIsNotNone(a1)
        b1 = self.h5file.create_array(root, "b1", b)
        self.assertIsNotNone(b1)
        expr = tb.Expr("a | b")
        r1 = expr.eval()
        r2 = eval("a | b")
        if common.verbose:
            print("Computed expression:", repr(r1), r1.dtype)
            print("Should look like:", repr(r2), r2.dtype)
        self.assertTrue(
            common.areArraysEqual(r1, r2),
            "Evaluate is returning a wrong value.",
        )

    def test01_shortint(self):
        """Checking int8,uint8,int16,uint16 and int32 in expression."""

        for dtype in "int8", "uint8", "int16", "uint16", "int32":
            if common.verbose:
                print("Checking type:", dtype)
            # Build arrays with different shapes as inputs
            a = np.array([1, 2, 3], dtype)
            b = np.array([3, 4, 5], dtype)
            root = self.h5file.root
            a1 = self.h5file.create_array(root, "a1", a)
            b1 = self.h5file.create_array(root, "b1", b)
            two = np.int32(2)
            self.assertIsInstance(two, np.integer)
            expr = tb.Expr("two * a1-b1")
            r1 = expr.eval()
            a = np.array([1, 2, 3], "int32")
            b = np.array([3, 4, 5], "int32")
            r2 = eval("two * a-b")
            if common.verbose:
                print("Computed expression:", repr(r1), r1.dtype)
                print("Should look like:", repr(r2), r2.dtype)
            self.assertEqual(r1.dtype, r2.dtype)
            self.assertTrue(
                common.areArraysEqual(r1, r2),
                "Evaluate is returning a wrong value.",
            )
            # Remove created leaves
            a1.remove()
            b1.remove()

    def test02_longint(self):
        """Checking uint32 and int64 in expression."""

        for dtype in "uint32", "int64":
            if common.verbose:
                print("Checking type:", dtype)
            # Build arrays with different shapes as inputs
            a = np.array([1, 2, 3], dtype)
            b = np.array([3, 4, 5], dtype)
            root = self.h5file.root
            a1 = self.h5file.create_array(root, "a1", a)
            b1 = self.h5file.create_array(root, "b1", b)
            expr = tb.Expr("2 * a1-b1")
            r1 = expr.eval()
            a = np.array([1, 2, 3], "int64")
            b = np.array([3, 4, 5], "int64")
            r2 = eval("2 * a-b")
            if common.verbose:
                print("Computed expression:", repr(r1), r1.dtype)
                print("Should look like:", repr(r2), r2.dtype)
            self.assertEqual(r1.dtype, r2.dtype)
            self.assertTrue(
                common.areArraysEqual(r1, r2),
                "Evaluate is returning a wrong value.",
            )
            # Remove created leaves
            a1.remove()
            b1.remove()

    def test03_float(self):
        """Checking float32 and float64 in expression."""

        for dtype in "float32", "float64":
            if common.verbose:
                print("Checking type:", dtype)
            # Build arrays with different shapes as inputs
            a = np.array([1, 2, 3], dtype)
            b = np.array([3, 4, 5], dtype)
            root = self.h5file.root
            a1 = self.h5file.create_array(root, "a1", a)
            b1 = self.h5file.create_array(root, "b1", b)
            expr = tb.Expr("2 * a1-b1")
            r1 = expr.eval()
            a = np.array([1, 2, 3], dtype)
            b = np.array([3, 4, 5], dtype)
            r2 = eval("2 * a-b")
            if common.verbose:
                print("Computed expression:", repr(r1), r1.dtype)
                print("Should look like:", repr(r2), r2.dtype)
            self.assertEqual(r1.dtype, r2.dtype)
            self.assertTrue(
                common.areArraysEqual(r1, r2),
                "Evaluate is returning a wrong value.",
            )
            # Remove created leaves
            a1.remove()
            b1.remove()

    def test04_complex(self):
        """Checking complex64 and complex128 in expression."""

        for dtype in "complex64", "complex128":
            if common.verbose:
                print("Checking type:", dtype)
            # Build arrays with different shapes as inputs
            a = np.array([1, 2j, 3 + 2j], dtype)
            b = np.array([3, 4j, 5 + 1j], dtype)
            root = self.h5file.root
            a1 = self.h5file.create_array(root, "a1", a)
            b1 = self.h5file.create_array(root, "b1", b)
            expr = tb.Expr("2 * a1-b1")
            r1 = expr.eval()
            a = np.array([1, 2j, 3 + 2j], "complex128")
            b = np.array([3, 4j, 5 + 1j], "complex128")
            r2 = eval("2 * a-b")
            if common.verbose:
                print("Computed expression:", repr(r1), r1.dtype)
                print("Should look like:", repr(r2), r2.dtype)
            self.assertEqual(r1.dtype, r2.dtype)
            self.assertTrue(
                common.areArraysEqual(r1, r2),
                "Evaluate is returning a wrong value.",
            )
            # Remove created leaves
            a1.remove()
            b1.remove()

    def test05_string(self):
        """Checking strings in expression."""

        # Build arrays with different shapes as inputs
        a = np.array(["a", "bd", "cd"], "S")
        b = np.array(["a", "bdcd", "ccdc"], "S")
        root = self.h5file.root
        a1 = self.h5file.create_array(root, "a1", a)
        self.assertIsNotNone(a1)
        b1 = self.h5file.create_array(root, "b1", b)
        self.assertIsNotNone(b1)
        expr = tb.Expr("(a1 > b'a') | ( b1 > b'b')")
        r1 = expr.eval()
        r2 = eval("(a > b'a') | ( b > b'b')")
        if common.verbose:
            print("Computed expression:", repr(r1), r1.dtype)
            print("Should look like:", repr(r2), r2.dtype)
        self.assertTrue(
            common.areArraysEqual(r1, r2),
            "Evaluate is returning a wrong value.",
        )


# Test for different functions
class FunctionsTestCase(common.TempFileMixin, common.PyTablesTestCase):

    def test00_simple(self):
        """Checking some math functions in expression."""

        # Build arrays with different shapes as inputs
        a = np.array([0.1, 0.2, 0.3])
        b = np.array([0.3, 0.4, 0.5])
        root = self.h5file.root
        a1 = self.h5file.create_array(root, "a1", a)
        self.assertIsNotNone(a1)
        b1 = self.h5file.create_array(root, "b1", b)
        self.assertIsNotNone(b1)
        # The expression
        expr = tb.Expr("sin(a1) * sqrt(b1)")
        r1 = expr.eval()
        r2 = np.sin(a) * np.sqrt(b)
        if common.verbose:
            print("Computed expression:", repr(r1), r1.dtype)
            print("Should look like:", repr(r2), r2.dtype)
        npt.assert_array_almost_equal_nulp(r1, r2)


# Test for EArrays with maindim != 0
class MaindimTestCase(common.TempFileMixin, common.PyTablesTestCase):

    def test00_simple(self):
        """Checking other dimensions than 0 as main dimension."""

        shape = list(self.shape)
        # Build input arrays
        a = np.arange(np.prod(shape), dtype="i4").reshape(shape)
        b = a.copy()
        c = a.copy()
        root = self.h5file.root
        shape[self.maindim] = 0
        a1 = self.h5file.create_earray(
            root, "a1", atom=tb.Int32Col(), shape=shape
        )
        b1 = self.h5file.create_earray(
            root, "b1", atom=tb.Int32Col(), shape=shape
        )
        c1 = self.h5file.create_earray(
            root, "c1", atom=tb.Int32Col(), shape=shape
        )
        a1.append(a)
        b1.append(b)
        c1.append(c)
        # The expression
        expr = tb.Expr("2 * a1 + b1-c1")
        r1 = expr.eval()
        r2 = eval("2 * a + b-c")
        if common.verbose:
            print("Tested shape:", shape)
            print("Computed expression:", repr(r1), r1.dtype)
            print("Should look like:", repr(r2), r2.dtype)
        self.assertTrue(
            common.areArraysEqual(r1, r2),
            "Evaluate is returning a wrong value.",
        )

    def test01_out(self):
        """Checking other dimensions than 0 as main dimension (out)"""

        shape = list(self.shape)
        # Build input arrays
        a = np.arange(np.prod(shape), dtype="i4").reshape(shape)
        b = a.copy()
        c = a.copy()
        root = self.h5file.root
        shape[self.maindim] = 0
        a1 = self.h5file.create_earray(
            root, "a1", atom=tb.Int32Col(), shape=shape
        )
        b1 = self.h5file.create_earray(
            root, "b1", atom=tb.Int32Col(), shape=shape
        )
        c1 = self.h5file.create_earray(
            root, "c1", atom=tb.Int32Col(), shape=shape
        )
        r1 = self.h5file.create_earray(
            root, "r1", atom=tb.Int32Col(), shape=shape
        )
        a1.append(a)
        b1.append(b)
        c1.append(c)
        r1.append(c)
        # The expression
        expr = tb.Expr("2 * a1 + b1-c1")
        expr.set_output(r1)
        expr.eval()
        r2 = eval("2 * a + b-c")
        if common.verbose:
            print("Tested shape:", shape)
            print("Computed expression:", repr(r1[:]), r1.dtype)
            print("Should look like:", repr(r2), r2.dtype)
        self.assertTrue(
            common.areArraysEqual(r1[:], r2),
            "Evaluate is returning a wrong value.",
        )

    def test02_diff_in_maindims(self):
        """Checking different main dimensions in inputs."""

        shape = list(self.shape)
        # Build input arrays
        a = np.arange(np.prod(shape), dtype="i4").reshape(shape)
        b = a.copy()
        c = a.copy()
        root = self.h5file.root
        shape2 = shape[:]
        shape[self.maindim] = 0
        shape2[0] = 0
        a1 = self.h5file.create_earray(
            root, "a1", atom=tb.Int32Col(), shape=shape
        )
        self.assertEqual(a1.maindim, self.maindim)
        b1 = self.h5file.create_earray(
            root, "b1", atom=tb.Int32Col(), shape=shape2
        )
        self.assertEqual(b1.maindim, 0)
        c1 = self.h5file.create_earray(
            root, "c1", atom=tb.Int32Col(), shape=shape
        )
        r1 = self.h5file.create_earray(
            root, "r1", atom=tb.Int32Col(), shape=shape
        )
        a1.append(a)
        b1.append(b)
        c1.append(c)
        r1.append(c)
        # The expression
        expr = tb.Expr("2 * a1 + b1-c1")
        r1 = expr.eval()
        r2 = eval("2 * a + b-c")
        if common.verbose:
            print("Tested shape:", shape)
            print("Computed expression:", repr(r1), r1.dtype)
            print("Should look like:", repr(r2), r2.dtype)
        self.assertTrue(
            common.areArraysEqual(r1, r2),
            "Evaluate is returning a wrong value.",
        )

    def test03_diff_in_out_maindims(self):
        """Checking different maindims in inputs and output."""

        shape = list(self.shape)
        # Build input arrays
        a = np.arange(np.prod(shape), dtype="i4").reshape(shape)
        b = a.copy()
        c = a.copy()
        root = self.h5file.root
        shape2 = shape[:]
        shape[self.maindim] = 0
        shape2[0] = 0
        a1 = self.h5file.create_earray(
            root, "a1", atom=tb.Int32Col(), shape=shape
        )
        self.assertEqual(a1.maindim, self.maindim)
        b1 = self.h5file.create_earray(
            root, "b1", atom=tb.Int32Col(), shape=shape
        )
        c1 = self.h5file.create_earray(
            root, "c1", atom=tb.Int32Col(), shape=shape
        )
        r1 = self.h5file.create_earray(
            root, "r1", atom=tb.Int32Col(), shape=shape2
        )
        self.assertEqual(r1.maindim, 0)
        a1.append(a)
        b1.append(b)
        c1.append(c)
        r1.append(c)
        # The expression
        expr = tb.Expr("2 * a1 + b1-c1")
        expr.set_output(r1)
        expr.eval()
        r2 = eval("2 * a + b-c")
        if common.verbose:
            print("Tested shape:", shape)
            print("Computed expression:", repr(r1[:]), r1.dtype)
            print("Should look like:", repr(r2), r2.dtype)
        self.assertTrue(
            common.areArraysEqual(r1[:], r2),
            "Evaluate is returning a wrong value.",
        )

    def test04_diff_in_out_maindims_lengths(self):
        """Checking different maindims and lengths in inputs and output."""

        shape = list(self.shape)
        # Build input arrays
        a = np.arange(np.prod(shape), dtype="i4").reshape(shape)
        b = a.copy()
        c = a.copy()
        root = self.h5file.root
        shape2 = shape[:]
        shape[self.maindim] = 0
        shape2[0] = 0
        a1 = self.h5file.create_earray(
            root, "a1", atom=tb.Int32Col(), shape=shape
        )
        self.assertEqual(a1.maindim, self.maindim)
        b1 = self.h5file.create_earray(
            root, "b1", atom=tb.Int32Col(), shape=shape
        )
        c1 = self.h5file.create_earray(
            root, "c1", atom=tb.Int32Col(), shape=shape
        )
        r1 = self.h5file.create_earray(
            root, "r1", atom=tb.Int32Col(), shape=shape2
        )
        self.assertEqual(r1.maindim, 0)
        a1.append(a)
        a1.append(a)
        b1.append(b)
        b1.append(b)
        c1.append(c)
        c1.append(c)
        r1.append(c)  # just once so that output is smaller
        # The expression
        expr = tb.Expr("2 * a1 + b1-c1")
        expr.set_output(r1)
        # This should raise an error
        self.assertRaises(ValueError, expr.eval)


class Maindim0(MaindimTestCase):
    maindim = 1
    shape = (1, 2)


class Maindim1(MaindimTestCase):
    maindim = 1
    shape = (2, 3)


class Maindim2(MaindimTestCase):
    maindim = 1
    shape = (2, 3, 4)


class Maindim3(MaindimTestCase):
    maindim = 2
    shape = (2, 3, 4)


# Test `append` mode flag in `set_output()`
class AppendModeTestCase(common.TempFileMixin, common.PyTablesTestCase):

    def test01_append(self):
        """Checking append mode in `set_output()`"""

        shape = [3, 2]
        # Build input arrays
        a = np.arange(np.prod(shape), dtype="i4").reshape(shape)
        b = a.copy()
        c = a.copy()
        shape[1] = 0
        root = self.h5file.root
        a1 = self.h5file.create_earray(
            root, "a1", atom=tb.Int32Col(), shape=shape
        )
        b1 = self.h5file.create_earray(
            root, "b1", atom=tb.Int32Col(), shape=shape
        )
        c1 = self.h5file.create_earray(
            root, "c1", atom=tb.Int32Col(), shape=shape
        )
        r1 = self.h5file.create_earray(
            root, "r1", atom=tb.Int32Col(), shape=shape
        )
        a1.append(a)
        b1.append(b)
        c1.append(c)
        if not self.append:
            r1.append(c)
        # The expression
        expr = tb.Expr("2 * a1 + b1-c1")
        expr.set_output(r1, append_mode=self.append)
        expr.eval()
        r2 = eval("2 * a + b-c")
        if common.verbose:
            print("Tested shape:", shape)
            print("Computed expression:", repr(r1[:]), r1.dtype)
            print("Should look like:", repr(r2), r2.dtype)
        self.assertTrue(
            common.areArraysEqual(r1[:], r2),
            "Evaluate is returning a wrong value.",
        )


class AppendModeTrue(AppendModeTestCase):
    append = True


class AppendModeFalse(AppendModeTestCase):
    append = False


# Test for `__iter__()` iterator
class iterTestCase(common.TempFileMixin, common.PyTablesTestCase):

    def setUp(self):
        super().setUp()
        shape = list(self.shape)
        # Build input arrays
        a = np.arange(np.prod(shape), dtype="i4").reshape(shape)
        b = a.copy()
        c = a.copy()
        self.npvars = {"a": a, "b": b, "c": c}
        shape[self.maindim] = 0
        root = self.h5file.root
        a1 = self.h5file.create_earray(
            root, "a1", atom=tb.Int32Col(), shape=shape
        )
        b1 = self.h5file.create_earray(
            root, "b1", atom=tb.Int32Col(), shape=shape
        )
        c1 = self.h5file.create_earray(
            root, "c1", atom=tb.Int32Col(), shape=shape
        )
        a1.append(a)
        b1.append(b)
        c1.append(c)
        self.vars = {"a": a1, "b": b1, "c": c1}
        # The expression
        self.sexpr = "2 * a + b-c"

    def test00_iter(self):
        """Checking the __iter__ iterator."""

        expr = tb.Expr(self.sexpr, self.vars)
        r1 = np.array([row for row in expr])
        r2 = eval(self.sexpr, self.npvars)
        if common.verbose:
            print("Tested shape, maindim:", self.shape, self.maindim)
            print("Computed expression:", repr(r1[:]), r1.dtype)
            print("Should look like:", repr(r2), r2.dtype)
        self.assertTrue(
            common.areArraysEqual(r1[:], r2),
            "Evaluate is returning a wrong value.",
        )

    def test01a_sss(self):
        """Checking the __iter__ iterator (with ranges, I)"""

        start, stop, step = self.range_[0], None, None
        expr = tb.Expr(self.sexpr, self.vars)
        expr.set_inputs_range(start, stop, step)
        r1 = np.array([row for row in expr])
        npvars = get_sliced_vars2(
            self.npvars, start, stop, step, self.shape, self.maindim
        )
        r2 = eval(self.sexpr, npvars)
        if common.verbose:
            print("Tested shape, maindim:", self.shape, self.maindim)
            print("Computed expression:", repr(r1[:]), r1.dtype)
            print("Should look like:", repr(r2), r2.dtype)
        self.assertTrue(
            common.areArraysEqual(r1[:], r2),
            "Evaluate is returning a wrong value.",
        )

    def test01b_sss(self):
        """Checking the __iter__ iterator (with ranges, II)"""

        start, stop, step = self.range_[0], self.range_[2], None
        expr = tb.Expr(self.sexpr, self.vars)
        expr.set_inputs_range(start, stop, step)
        r1 = np.array([row for row in expr])
        npvars = get_sliced_vars2(
            self.npvars, start, stop, step, self.shape, self.maindim
        )
        r2 = eval(self.sexpr, npvars)
        if common.verbose:
            print("Tested shape, maindim:", self.shape, self.maindim)
            print("Computed expression:", repr(r1[:]), r1.dtype)
            print("Should look like:", repr(r2), r2.dtype)
        self.assertTrue(
            common.areArraysEqual(r1[:], r2),
            "Evaluate is returning a wrong value.",
        )

    def test01c_sss(self):
        """Checking the __iter__ iterator (with ranges, III)"""

        start, stop, step = self.range_
        expr = tb.Expr(self.sexpr, self.vars)
        expr.set_inputs_range(start, stop, step)
        r1 = np.array([row for row in expr])
        npvars = get_sliced_vars2(
            self.npvars, start, stop, step, self.shape, self.maindim
        )
        r2 = eval(self.sexpr, npvars)
        if common.verbose:
            print("Tested shape, maindim:", self.shape, self.maindim)
            print("Computed expression:", repr(r1[:]), r1.dtype)
            print("Should look like:", repr(r2), r2.dtype)
        self.assertTrue(
            common.areArraysEqual(r1[:], r2),
            "Evaluate is returning a wrong value.",
        )


class iter0(iterTestCase):
    maindim = 0
    shape = (0,)
    range_ = (1, 2, 1)


class iter1(iterTestCase):
    maindim = 0
    shape = (3,)
    range_ = (1, 2, 1)


class iter2(iterTestCase):
    maindim = 0
    shape = (3, 2)
    range_ = (0, 3, 2)


class iter3(iterTestCase):
    maindim = 1
    shape = (3, 2)
    range_ = (0, 3, 2)


class iter4(iterTestCase):
    maindim = 2
    shape = (3, 2, 1)
    range_ = (1, 3, 2)


class iter5(iterTestCase):
    maindim = 2
    shape = (1, 2, 5)
    range_ = (0, 4, 2)


# Test for set_output_range
class setOutputRangeTestCase(common.TempFileMixin, common.PyTablesTestCase):

    def test00_simple(self):
        """Checking the range selection for output."""

        shape = list(self.shape)
        start, stop, step = self.range_
        # Build input arrays
        a = np.arange(np.prod(shape), dtype="i4").reshape(shape)
        b = a.copy()
        r = a.copy()
        root = self.h5file.root
        a1 = self.h5file.create_array(root, "a1", a)
        self.assertIsNotNone(a1)
        b1 = self.h5file.create_array(root, "b1", b)
        self.assertIsNotNone(b1)
        r1 = self.h5file.create_array(root, "r1", r)
        # The expression
        expr = tb.Expr("a1-b1-1")
        expr.set_output(r1)
        expr.set_output_range(start, stop, step)
        expr.eval()
        r2 = eval("a-b-1")
        r[start:stop:step] = r2[: len(range(start, stop, step))]
        if common.verbose:
            print("Tested shape:", shape)
            print("Computed expression:", repr(r1[:]), r1.dtype)
            print("Should look like:", repr(r), r.dtype)
        self.assertTrue(
            common.areArraysEqual(r1[:], r),
            "Evaluate is returning a wrong value.",
        )

    def test01_maindim(self):
        """Checking the range selection for output (maindim > 0)"""

        shape = list(self.shape)
        start, stop, step = self.range_
        # Build input arrays
        a = np.arange(np.prod(shape), dtype="i4").reshape(shape)
        b = a.copy()
        r = a.copy()
        shape[self.maindim] = 0
        root = self.h5file.root
        a1 = self.h5file.create_earray(
            root, "a1", atom=tb.Int32Col(), shape=shape
        )
        b1 = self.h5file.create_earray(
            root, "b1", atom=tb.Int32Col(), shape=shape
        )
        r1 = self.h5file.create_earray(
            root, "r1", atom=tb.Int32Col(), shape=shape
        )
        a1.append(a)
        b1.append(b)
        r1.append(r)
        # The expression
        expr = tb.Expr("a1-b1-1")
        expr.set_output(r1)
        expr.set_output_range(start, stop, step)
        expr.eval()
        r2 = eval("a-b-1")
        lsl = tuple([slice(None)] * self.maindim)
        # print "lsl-->", lsl + (slice(start,stop,step),)
        lrange = len(range(start, stop, step))
        r.__setitem__(
            lsl + (slice(start, stop, step),),
            r2.__getitem__(lsl + (slice(0, lrange),)),
        )
        if common.verbose:
            print("Tested shape:", shape)
            print("Computed expression:", repr(r1[:]), r1.dtype)
            print("Should look like:", repr(r), r.dtype)
        self.assertTrue(
            common.areArraysEqual(r1[:], r),
            "Evaluate is returning a wrong value.",
        )


class setOutputRange0(setOutputRangeTestCase):
    maindim = 0
    shape = (10,)
    range_ = (0, 1, 2)


class setOutputRange1(setOutputRangeTestCase):
    maindim = 0
    shape = (10,)
    range_ = (0, 10, 2)


class setOutputRange2(setOutputRangeTestCase):
    maindim = 0
    shape = (10,)
    range_ = (1, 10, 2)


class setOutputRange3(setOutputRangeTestCase):
    maindim = 0
    shape = (10, 1)
    range_ = (1, 10, 3)


class setOutputRange4(setOutputRangeTestCase):
    maindim = 0
    shape = (10, 2)
    range_ = (1, 10, 3)


class setOutputRange5(setOutputRangeTestCase):
    maindim = 0
    shape = (5, 3, 1)
    range_ = (1, 5, 1)


class setOutputRange6(setOutputRangeTestCase):
    maindim = 1
    shape = (2, 5)
    range_ = (1, 3, 2)


class setOutputRange7(setOutputRangeTestCase):
    maindim = 1
    shape = (2, 5, 1)
    range_ = (1, 3, 2)


class setOutputRange8(setOutputRangeTestCase):
    maindim = 2
    shape = (1, 3, 5)
    range_ = (1, 5, 2)


class setOutputRange9(setOutputRangeTestCase):
    maindim = 3
    shape = (1, 3, 4, 5)
    range_ = (1, 5, 3)


# Test for very large inputs
class VeryLargeInputsTestCase(common.TempFileMixin, common.PyTablesTestCase):

    def test00_simple(self):
        """Checking very large inputs."""

        shape = self.shape
        # Use filters so as to not use too much space
        if tb.which_lib_version("blosc") is not None:
            filters = tb.Filters(complevel=1, complib="blosc", shuffle=False)
        elif tb.which_lib_version("lzo") is not None:
            filters = tb.Filters(complevel=1, complib="lzo", shuffle=False)
        else:
            filters = tb.Filters(complevel=1, shuffle=False)
        # Build input arrays
        root = self.h5file.root
        a = self.h5file.create_carray(
            root,
            "a",
            atom=tb.Float64Atom(dflt=3),
            shape=shape,
            filters=filters,
        )
        self.assertIsNotNone(a)
        b = self.h5file.create_carray(
            root,
            "b",
            atom=tb.Float64Atom(dflt=2),
            shape=shape,
            filters=filters,
        )
        self.assertIsNotNone(b)
        r1 = self.h5file.create_carray(
            root,
            "r1",
            atom=tb.Float64Atom(dflt=3),
            shape=shape,
            filters=filters,
        )
        # The expression
        expr = tb.Expr("a * b-6")  # Should give 0
        expr.set_output(r1)
        expr.eval()
        r1 = r1[-10:]  # Get the last ten rows
        r2 = np.zeros(10, dtype="float64")
        if common.verbose:
            print("Tested shape:", shape)
            print("Ten last rows:", repr(r1), r1.dtype)
            print("Should look like:", repr(r2), r2.dtype)
        self.assertTrue(
            common.areArraysEqual(r1, r2),
            "Evaluate is returning a wrong value.",
        )

    def test01_iter(self):
        """Checking very large inputs (__iter__ version)"""

        shape = self.shape
        if shape[0] >= 2**24:
            # The iterator is much slower, so don't run it for
            # extremely large arrays.
            if common.verbose:
                print("Skipping this *very* long test")
            return
        # Use filters so as to not use too much space
        if tb.which_lib_version("lzo") is not None:
            filters = tb.Filters(complevel=1, complib="lzo", shuffle=False)
        else:
            filters = tb.Filters(complevel=1, shuffle=False)

        # Build input arrays
        root = self.h5file.root
        a = self.h5file.create_carray(
            root, "a", atom=tb.Int32Atom(dflt=1), shape=shape, filters=filters
        )
        self.assertIsNotNone(a)
        b = self.h5file.create_carray(
            root, "b", atom=tb.Int32Atom(dflt=2), shape=shape, filters=filters
        )
        self.assertIsNotNone(b)
        r1 = self.h5file.create_carray(
            root, "r1", atom=tb.Int32Atom(dflt=3), shape=shape, filters=filters
        )
        # The expression
        expr = tb.Expr("a-b + 1")
        r1 = sum(expr)  # Should give 0
        if common.verbose:
            print("Tested shape:", shape)
            print("Cummulated sum:", r1)
            print("Should look like:", 0)
        self.assertEqual(r1, 0, "Evaluate is returning a wrong value.")


# The next can go on regular tests, as it should be light enough
class VeryLargeInputs1(VeryLargeInputsTestCase):
    shape = (2**20,)  # larger than any internal I/O buffers


# The next is only meant for 'heavy' mode as it can take more than 1 minute
# on modern machines
class VeryLargeInputs2(VeryLargeInputsTestCase):
    shape = (2**32 + 1,)  # check that arrays > 32-bit are supported


def suite():
    """Return a test suite consisting of all the test cases in the module."""

    theSuite = common.unittest.TestSuite()
    niter = 1
    # common.heavy = 1  # uncomment this only for testing purposes

    for i in range(niter):
        theSuite.addTest(common.make_suite(ExprNumPy))
        theSuite.addTest(common.make_suite(ExprArray))
        theSuite.addTest(common.make_suite(ExprCArray))
        theSuite.addTest(common.make_suite(ExprEArray))
        theSuite.addTest(common.make_suite(ExprColumn))
        theSuite.addTest(common.make_suite(MixedContainers0))
        theSuite.addTest(common.make_suite(MixedContainers1))
        theSuite.addTest(common.make_suite(MixedContainers2))
        theSuite.addTest(common.make_suite(MixedContainers3))
        theSuite.addTest(common.make_suite(UnalignedObject))
        theSuite.addTest(common.make_suite(NonContiguousObject))
        theSuite.addTest(common.make_suite(ExprError))
        theSuite.addTest(common.make_suite(Broadcast0))
        theSuite.addTest(common.make_suite(Broadcast1))
        theSuite.addTest(common.make_suite(Broadcast2))
        theSuite.addTest(common.make_suite(Broadcast3))
        theSuite.addTest(common.make_suite(Broadcast4))
        theSuite.addTest(common.make_suite(Broadcast5))
        theSuite.addTest(common.make_suite(DiffLength0))
        theSuite.addTest(common.make_suite(DiffLength1))
        theSuite.addTest(common.make_suite(DiffLength2))
        theSuite.addTest(common.make_suite(DiffLength3))
        theSuite.addTest(common.make_suite(DiffLength4))
        theSuite.addTest(common.make_suite(TypesTestCase))
        theSuite.addTest(common.make_suite(FunctionsTestCase))
        theSuite.addTest(common.make_suite(Maindim0))
        theSuite.addTest(common.make_suite(Maindim1))
        theSuite.addTest(common.make_suite(Maindim2))
        theSuite.addTest(common.make_suite(Maindim3))
        theSuite.addTest(common.make_suite(AppendModeTrue))
        theSuite.addTest(common.make_suite(AppendModeFalse))
        theSuite.addTest(common.make_suite(iter0))
        theSuite.addTest(common.make_suite(iter1))
        theSuite.addTest(common.make_suite(iter2))
        theSuite.addTest(common.make_suite(iter3))
        theSuite.addTest(common.make_suite(iter4))
        theSuite.addTest(common.make_suite(iter5))
        theSuite.addTest(common.make_suite(setOutputRange0))
        theSuite.addTest(common.make_suite(setOutputRange1))
        theSuite.addTest(common.make_suite(setOutputRange2))
        theSuite.addTest(common.make_suite(setOutputRange3))
        theSuite.addTest(common.make_suite(setOutputRange4))
        theSuite.addTest(common.make_suite(setOutputRange5))
        theSuite.addTest(common.make_suite(setOutputRange6))
        theSuite.addTest(common.make_suite(setOutputRange7))
        theSuite.addTest(common.make_suite(setOutputRange8))
        theSuite.addTest(common.make_suite(setOutputRange9))
        theSuite.addTest(common.make_suite(VeryLargeInputs1))
        if common.heavy:
            theSuite.addTest(common.make_suite(VeryLargeInputs2))
    return theSuite


if __name__ == "__main__":
    import sys

    common.parse_argv(sys.argv)
    common.print_versions()
    common.unittest.main(defaultTest="suite")
