"""Here is defined the Expr class."""

import sys
import warnings

import numexpr as ne
import numpy as np
import tables as tb

from .exceptions import PerformanceWarning
from .parameters import IO_BUFFER_SIZE, BUFFER_TIMES


class Expr:
    """A class for evaluating expressions with arbitrary array-like objects.

    Expr is a class for evaluating expressions containing array-like objects.
    With it, you can evaluate expressions (like "3 * a + 4 * b") that
    operate on arbitrary large arrays while optimizing the resources
    required to perform them (basically main memory and CPU cache memory).
    It is similar to the Numexpr package (see :ref:`[NUMEXPR] <NUMEXPR>`),
    but in addition to NumPy objects, it also accepts disk-based homogeneous
    arrays, like the Array, CArray, EArray and Column PyTables objects.

    .. warning::

        Expr class only offers a subset of the Numexpr features due to the
        complexity of implement some of them when dealing with huge amount of
        data.

    All the internal computations are performed via the Numexpr package,
    so all the broadcast and upcasting rules of Numexpr applies here too.
    These rules are very similar to the NumPy ones, but with some exceptions
    due to the particularities of having to deal with potentially very large
    disk-based arrays.  Be sure to read the documentation of the Expr
    constructor and methods as well as that of Numexpr, if you want to fully
    grasp these particularities.


    Parameters
    ----------
    expr : str
        This specifies the expression to be evaluated, such as "2 * a + 3 * b".
    uservars : dict
        This can be used to define the variable names appearing in *expr*.
        This mapping should consist of identifier-like strings pointing to any
        `Array`, `CArray`, `EArray`, `Column` or NumPy ndarray instances (or
        even others which will tried to be converted to ndarrays).  When
        `uservars` is not provided or `None`, the current local and global
        namespace is sought instead of `uservars`.  It is also possible to pass
        just some of the variables in expression via the `uservars` mapping,
        and the rest will be retrieved from the current local and global
        namespaces.
    kwargs : dict
        This is meant to pass additional parameters to the Numexpr kernel.
        This is basically the same as the kwargs argument in
        Numexpr.evaluate(), and is mainly meant for advanced use.

    Examples
    --------
    The following shows an example of using Expr::

        >>> f = tb.open_file('/tmp/test_expr.h5', 'w')
        >>> a = f.create_array('/', 'a', np.array([1,2,3]))
        >>> b = f.create_array('/', 'b', np.array([3,4,5]))
        >>> c = np.array([4,5,6])
        >>> expr = tb.Expr("2 * a + b * c")   # initialize the expression
        >>> expr.eval()                 # evaluate it
        array([14, 24, 36], dtype=int64)
        >>> sum(expr)                   # use as an iterator
        74

    where you can see that you can mix different containers in
    the expression (whenever shapes are consistent).

    You can also work with multidimensional arrays::

        >>> a2 = f.create_array('/', 'a2', np.array([[1,2],[3,4]]))
        >>> b2 = f.create_array('/', 'b2', np.array([[3,4],[5,6]]))
        >>> c2 = np.array([4,5])           # This will be broadcasted
        >>> expr = tb.Expr("2 * a2 + b2-c2")
        >>> expr.eval()
        array([[1, 3],
               [7, 9]], dtype=int64)
        >>> sum(expr)
        array([ 8, 12], dtype=int64)
        >>> f.close()

    .. rubric:: Expr attributes

    .. attribute:: append_mode

        The append mode for user-provided output containers.

    .. attribute:: maindim

        Common main dimension for inputs in expression.

    .. attribute:: names

        The names of variables in expression (list).

    .. attribute:: out

        The user-provided container (if any) for the expression outcome.

    .. attribute:: o_start

        The start range selection for the user-provided output.

    .. attribute:: o_stop

        The stop range selection for the user-provided output.

    .. attribute:: o_step

        The step range selection for the user-provided output.

    .. attribute:: shape

        Common shape for the arrays in expression.

    .. attribute:: values

        The values of variables in expression (list).

    """

    _exprvars_cache = {}
    """Cache of variables participating in expressions.

    .. versionadded:: 3.0

    """

    def __init__(self, expr, uservars=None, **kwargs):

        self.append_mode = False
        """The append mode for user-provided output containers."""
        self.maindim = 0
        """Common main dimension for inputs in expression."""
        self.names = []
        """The names of variables in expression (list)."""
        self.out = None
        """The user-provided container (if any) for the expression outcome."""
        self.o_start = None
        """The start range selection for the user-provided output."""
        self.o_stop = None
        """The stop range selection for the user-provided output."""
        self.o_step = None
        """The step range selection for the user-provided output."""
        self.shape = None
        """Common shape for the arrays in expression."""
        self.start, self.stop, self.step = (None,) * 3
        self.start = None
        """The start range selection for the input."""
        self.stop = None
        """The stop range selection for the input."""
        self.step = None
        """The step range selection for the input."""
        self.values = []
        """The values of variables in expression (list)."""

        self._compiled_expr = None
        """The compiled expression."""
        self._single_row_out = None
        """A sample of the output with just a single row."""

        # First, get the signature for the arrays in expression
        vars_ = self._required_expr_vars(expr, uservars)
        context = ne.necompiler.getContext(kwargs)
        self.names, _ = ne.necompiler.getExprNames(expr, context)

        # Raise a ValueError in case we have unsupported objects
        for name, var in vars_.items():
            if type(var) in (int, float, str):
                continue
            if not isinstance(var, (tb.Leaf, tb.Column)):
                if hasattr(var, "dtype"):
                    # Quacks like a NumPy object
                    continue
                raise TypeError("Unsupported variable type: %r" % var)
            objname = var.__class__.__name__
            if objname not in ("Array", "CArray", "EArray", "Column"):
                raise TypeError("Unsupported variable type: %r" % var)

        # NumPy arrays to be copied? (we don't need to worry about
        # PyTables objects, as the reads always return contiguous and
        # aligned objects, or at least I think so).
        for name, var in vars_.items():
            if isinstance(var, np.ndarray):
                # See numexpr.necompiler.evaluate for a rational
                # of the code below
                if not var.flags.aligned:
                    if var.ndim != 1:
                        # Do a copy of this variable
                        var = var.copy()
                        # Update the vars_ dictionary
                        vars_[name] = var

        # Get the variables and types
        values = self.values
        types_ = []
        for name in self.names:
            value = vars_[name]
            if hasattr(value, 'atom'):
                types_.append(value.atom)
            elif hasattr(value, 'dtype'):
                types_.append(value)
            else:
                # try to convert into a NumPy array
                value = np.array(value)
                types_.append(value)
            values.append(value)

        # Create a signature for the expression
        signature = [(name, ne.necompiler.getType(type_))
                     for (name, type_) in zip(self.names, types_)]

        # Compile the expression
        self._compiled_expr = ne.necompiler.NumExpr(expr, signature, **kwargs)

        # Guess the shape for the outcome and the maindim of inputs
        self.shape, self.maindim = self._guess_shape()

    # The next method is similar to their counterpart in `Table`, but
    # adapted to the `Expr` own requirements.
    def _required_expr_vars(self, expression, uservars, depth=2):
        """Get the variables required by the `expression`.

        A new dictionary defining the variables used in the `expression`
        is returned.  Required variables are first looked up in the
        `uservars` mapping, then in the set of top-level columns of the
        table.  Unknown variables cause a `NameError` to be raised.

        When `uservars` is `None`, the local and global namespace where
        the API callable which uses this method is called is sought
        instead.  To disable this mechanism, just specify a mapping as
        `uservars`.

        Nested columns and variables with an ``uint64`` type are not
        allowed (`TypeError` and `NotImplementedError` are raised,
        respectively).

        `depth` specifies the depth of the frame in order to reach local
        or global variables.

        """

        # Get the names of variables used in the expression.
        exprvars_cache = self._exprvars_cache
        if expression not in exprvars_cache:
            # Protection against growing the cache too much
            if len(exprvars_cache) > 256:
                # Remove 10 (arbitrary) elements from the cache
                for k in list(exprvars_cache)[:10]:
                    del exprvars_cache[k]
            cexpr = compile(expression, '<string>', 'eval')
            exprvars = [var for var in cexpr.co_names
                        if var not in ['None', 'False', 'True']
                        and var not in ne.expressions.functions]
            exprvars_cache[expression] = exprvars
        else:
            exprvars = exprvars_cache[expression]

        # Get the local and global variable mappings of the user frame
        # if no mapping has been explicitly given for user variables.
        user_locals, user_globals = {}, {}
        if uservars is None:
            user_frame = sys._getframe(depth)
            user_locals = user_frame.f_locals
            user_globals = user_frame.f_globals

        # Look for the required variables first among the ones
        # explicitly provided by the user.
        reqvars = {}
        for var in exprvars:
            # Get the value.
            if uservars is not None and var in uservars:
                val = uservars[var]
            elif uservars is None and var in user_locals:
                val = user_locals[var]
            elif uservars is None and var in user_globals:
                val = user_globals[var]
            else:
                raise NameError("name ``%s`` is not defined" % var)

            # Check the value.
            if hasattr(val, 'dtype') and val.dtype.str[1:] == 'u8':
                raise NotImplementedError(
                    "variable ``%s`` refers to "
                    "a 64-bit unsigned integer object, that is "
                    "not yet supported in expressions, sorry; " % var)
            elif hasattr(val, '_v_colpathnames'):  # nested column
                # This branch is never reached because the compile step
                # above already raise a ``TypeError`` for nested
                # columns, but that could change in the future.  So it
                # is best to let this here.
                raise TypeError(
                    "variable ``%s`` refers to a nested column, "
                    "not allowed in expressions" % var)
            reqvars[var] = val
        return reqvars

    def set_inputs_range(self, start=None, stop=None, step=None):
        """Define a range for all inputs in expression.

        The computation will only take place for the range defined by
        the start, stop and step parameters in the main dimension of
        inputs (or the leading one, if the object lacks the concept of
        main dimension, like a NumPy container).  If not a common main
        dimension exists for all inputs, the leading dimension will be
        used instead.

        """

        self.start = start
        self.stop = stop
        self.step = step

    def set_output(self, out, append_mode=False):
        """Set out as container for output as well as the append_mode.

        The out must be a container that is meant to keep the outcome of
        the expression.  It should be an homogeneous type container and
        can typically be an Array, CArray, EArray, Column or a NumPy ndarray.

        The append_mode specifies the way of which the output is filled.
        If true, the rows of the outcome are *appended* to the out container.
        Of course, for doing this it is necessary that out would have an
        append() method (like an EArray, for example).

        If append_mode is false, the output is set via the __setitem__()
        method (see the Expr.set_output_range() for info on how to select
        the rows to be updated).  If out is smaller than what is required
        by the expression, only the computations that are needed to fill
        up the container are carried out.  If it is larger, the excess
        elements are unaffected.

        """

        if not (hasattr(out, "shape") and hasattr(out, "__setitem__")):
            raise ValueError(
                "You need to pass a settable multidimensional container "
                "as output")
        self.out = out
        if append_mode and not hasattr(out, "append"):
            raise ValueError(
                "For activating the ``append`` mode, you need a container "
                "with an `append()` method (like the `EArray`)")
        self.append_mode = append_mode

    def set_output_range(self, start=None, stop=None, step=None):
        """Define a range for user-provided output object.

        The output object will only be modified in the range specified by the
        start, stop and step parameters in the main dimension of output (or the
        leading one, if the object does not have the concept of main dimension,
        like a NumPy container).

        """

        if self.out is None:
            raise IndexError(
                "You need to pass an output object to `setOut()` first")
        self.o_start = start
        self.o_stop = stop
        self.o_step = step

    # Although the next code is similar to the method in `Leaf`, it
    # allows the use of pure NumPy objects.
    def _calc_nrowsinbuf(self, object_):
        """Calculate the number of rows that will fit in a buffer."""

        # Compute the rowsize for the *leading* dimension
        shape_ = list(object_.shape)
        if shape_:
            shape_[0] = 1

        rowsize = np.prod(shape_) * object_.dtype.itemsize

        # Compute the nrowsinbuf
        # Multiplying the I/O buffer size by 4 gives optimal results
        # in my benchmarks with `tables.Expr` (see ``bench/poly.py``)
        buffersize = IO_BUFFER_SIZE * 4
        nrowsinbuf = buffersize // rowsize

        # Safeguard against row sizes being extremely large
        if nrowsinbuf == 0:
            nrowsinbuf = 1
            # If rowsize is too large, issue a Performance warning
            maxrowsize = BUFFER_TIMES * buffersize
            if rowsize > maxrowsize:
                warnings.warn("""\
The object ``%s`` is exceeding the maximum recommended rowsize (%d
bytes); be ready to see PyTables asking for *lots* of memory and
possibly slow I/O.  You may want to reduce the rowsize by trimming the
value of dimensions that are orthogonal (and preferably close) to the
*leading* dimension of this object."""
                              % (object, maxrowsize),
                              PerformanceWarning)

        return nrowsinbuf

    def _guess_shape(self):
        """Guess the shape of the output of the expression."""

        # First, compute the maximum dimension of inputs and maindim
        # (if it exists)
        maxndim = 0
        maindims = []
        for val in self.values:
            # Get the minimum of the lengths
            if len(val.shape) > maxndim:
                maxndim = len(val.shape)
            if hasattr(val, "maindim"):
                maindims.append(val.maindim)
        if maxndim == 0:
            self._single_row_out = out = self._compiled_expr(*self.values)
            return (), None
        if maindims and [maindims[0]] * len(maindims) == maindims:
            # If all maindims detected are the same, use this as maindim
            maindim = maindims[0]
        else:
            # If not, the main dimension will be the default one
            maindim = 0

        # The slices parameter for inputs
        slices = (slice(None),) * maindim + (0,)

        # Now, collect the values in first row of arrays with maximum dims
        vals = []
        lens = []
        for val in self.values:
            shape = val.shape
            # Warning: don't use len(val) below or it will raise an
            # `Overflow` error on 32-bit platforms for large enough arrays.
            if shape != () and shape[maindim] == 0:
                vals.append(val[:])
                lens.append(0)
            elif len(shape) < maxndim:
                vals.append(val)
            else:
                vals.append(val.__getitem__(slices))
                lens.append(shape[maindim])
        minlen = min(lens)
        self._single_row_out = out = self._compiled_expr(*vals)
        shape = list(out.shape)
        if minlen > 0:
            shape.insert(maindim, minlen)
        return shape, maindim

    def _get_info(self, shape, maindim, itermode=False):
        """Return various info needed for evaluating the computation loop."""

        # Compute the shape of the resulting container having
        # in account new possible values of start, stop and step in
        # the inputs range
        if maindim is not None:
            (start, stop, step) = slice(
                self.start, self.stop, self.step).indices(shape[maindim])
            shape[maindim] = min(
                shape[maindim], len(range(start, stop, step)))
            i_nrows = shape[maindim]
        else:
            start, stop, step = 0, 0, None
            i_nrows = 0

        if not itermode:
            # Create a container for output if not defined yet
            o_maindim = 0    # Default maindim
            if self.out is None:
                out = np.empty(shape, dtype=self._single_row_out.dtype)
                # Get the trivial values for start, stop and step
                if maindim is not None:
                    (o_start, o_stop, o_step) = (0, shape[maindim], 1)
                else:
                    (o_start, o_stop, o_step) = (0, 0, 1)
            else:
                out = self.out
                # Out container already provided.  Do some sanity checks.
                if hasattr(out, "maindim"):
                    o_maindim = out.maindim

                # Refine the shape of the resulting container having in
                # account new possible values of start, stop and step in
                # the output range
                o_shape = list(out.shape)
                s = slice(self.o_start, self.o_stop, self.o_step)
                o_start, o_stop, o_step = s.indices(o_shape[o_maindim])
                o_shape[o_maindim] = min(o_shape[o_maindim],
                                         len(range(o_start, o_stop, o_step)))

                # Check that the shape of output is consistent with inputs
                tr_oshape = list(o_shape)   # this implies a copy
                olen_ = tr_oshape.pop(o_maindim)
                tr_shape = list(shape)      # do a copy
                if maindim is not None:
                    len_ = tr_shape.pop(o_maindim)
                else:
                    len_ = 1
                if tr_oshape != tr_shape:
                    raise ValueError(
                        "Shape for out container does not match expression")
                # Force the input length to fit in `out`
                if not self.append_mode and olen_ < len_:
                    shape[o_maindim] = olen_
                    stop = start + olen_

        # Get the positions of inputs that should be sliced (the others
        # will be broadcasted)
        ndim = len(shape)
        slice_pos = [i for i, val in enumerate(self.values)
                     if len(val.shape) == ndim]

        # The size of the I/O buffer
        nrowsinbuf = 1
        for i, val in enumerate(self.values):
            # Skip scalar values in variables
            if i in slice_pos:
                nrows = self._calc_nrowsinbuf(val)
                if nrows > nrowsinbuf:
                    nrowsinbuf = nrows

        if not itermode:
            return (i_nrows, slice_pos, start, stop, step, nrowsinbuf,
                    out, o_maindim, o_start, o_stop, o_step)
        else:
            # For itermode, we don't need the out info
            return (i_nrows, slice_pos, start, stop, step, nrowsinbuf)

    def eval(self):
        """Evaluate the expression and return the outcome.

        Because of performance reasons, the computation order tries to go along
        the common main dimension of all inputs.  If not such a common main
        dimension is found, the iteration will go along the leading dimension
        instead.

        For non-consistent shapes in inputs (i.e. shapes having a different
        number of dimensions), the regular NumPy broadcast rules applies.
        There is one exception to this rule though: when the dimensions
        orthogonal to the main dimension of the expression are consistent, but
        the main dimension itself differs among the inputs, then the shortest
        one is chosen for doing the computations.  This is so because trying to
        expand very large on-disk arrays could be too expensive or simply not
        possible.

        Also, the regular Numexpr casting rules (which are similar to those of
        NumPy, although you should check the Numexpr manual for the exceptions)
        are applied to determine the output type.

        Finally, if the setOuput() method specifying a user container has
        already been called, the output is sent to this user-provided
        container.  If not, a fresh NumPy container is returned instead.

        .. warning::

            When dealing with large on-disk inputs, failing to specify an
            on-disk container may consume all your available memory.

        """

        values, shape, maindim = self.values, self.shape, self.maindim

        # Get different info we need for the main computation loop
        (i_nrows, slice_pos, start, stop, step, nrowsinbuf,
         out, o_maindim, o_start, o_stop, o_step) = \
            self._get_info(shape, maindim)

        if i_nrows == 0:
            # No elements to compute
            if start >= stop and self.start is not None:
                return out
            else:
                return self._single_row_out

        # Create a key that selects every element in inputs and output
        # (including the main dimension)
        i_slices = [slice(None)] * (maindim + 1)
        o_slices = [slice(None)] * (o_maindim + 1)

        # This is a hack to prevent doing unnecessary flavor conversions
        # while reading buffers
        for val in values:
            if hasattr(val, 'maindim'):
                val._v_convert = False

        # Start the computation itself
        for start2 in range(start, stop, step * nrowsinbuf):
            stop2 = start2 + step * nrowsinbuf
            if stop2 > stop:
                stop2 = stop
            # Set the proper slice for inputs
            i_slices[maindim] = slice(start2, stop2, step)
            # Get the input values
            vals = []
            for i, val in enumerate(values):
                if i in slice_pos:
                    vals.append(val.__getitem__(tuple(i_slices)))
                else:
                    # A read of values is not apparently needed, as PyTables
                    # leaves seems to work just fine inside Numexpr
                    vals.append(val)
            # Do the actual computation for this slice
            rout = self._compiled_expr(*vals)
            # Set the values into the out buffer
            if self.append_mode:
                out.append(rout)
            else:
                # Compute the slice to be filled in output
                start3 = o_start + (start2 - start) // step
                stop3 = start3 + nrowsinbuf * o_step
                if stop3 > o_stop:
                    stop3 = o_stop
                o_slices[o_maindim] = slice(start3, stop3, o_step)
                # Set the slice
                out[tuple(o_slices)] = rout

        # Activate the conversion again (default)
        for val in values:
            if hasattr(val, 'maindim'):
                val._v_convert = True

        return out

    def __iter__(self):
        """Iterate over the rows of the outcome of the expression.

        This iterator always returns rows as NumPy objects, so a possible out
        container specified in :meth:`Expr.set_output` method is ignored here.

        """

        values, shape, maindim = self.values, self.shape, self.maindim

        # Get different info we need for the main computation loop
        (i_nrows, slice_pos, start, stop, step, nrowsinbuf) = \
            self._get_info(shape, maindim, itermode=True)

        if i_nrows == 0:
            # No elements to compute
            return

        # Create a key that selects every element in inputs
        # (including the main dimension)
        i_slices = [slice(None)] * (maindim + 1)

        # This is a hack to prevent doing unnecessary flavor conversions
        # while reading buffers
        for val in values:
            if hasattr(val, 'maindim'):
                val._v_convert = False

        # Start the computation itself
        for start2 in range(start, stop, step * nrowsinbuf):
            stop2 = start2 + step * nrowsinbuf
            if stop2 > stop:
                stop2 = stop
            # Set the proper slice in the main dimension
            i_slices[maindim] = slice(start2, stop2, step)
            # Get the values for computing the buffer
            vals = []
            for i, val in enumerate(values):
                if i in slice_pos:
                    vals.append(val.__getitem__(tuple(i_slices)))
                else:
                    # A read of values is not apparently needed, as PyTables
                    # leaves seems to work just fine inside Numexpr
                    vals.append(val)
            # Do the actual computation
            rout = self._compiled_expr(*vals)
            # Return one row per call
            yield from rout

        # Activate the conversion again (default)
        for val in values:
            if hasattr(val, 'maindim'):
                val._v_convert = True


if __name__ == "__main__":

    # shape = (10000,10000)
    shape = (10, 10_000)

    f = tb.open_file("/tmp/expression.h5", "w")

    # Create some arrays
    a = f.create_carray(f.root, 'a', atom=tb.Float32Atom(dflt=1), shape=shape)
    b = f.create_carray(f.root, 'b', atom=tb.Float32Atom(dflt=2), shape=shape)
    c = f.create_carray(f.root, 'c', atom=tb.Float32Atom(dflt=3), shape=shape)
    out = f.create_carray(f.root, 'out', atom=tb.Float32Atom(dflt=3),
                          shape=shape)

    expr = Expr("a * b + c")
    expr.set_output(out)
    d = expr.eval()

    print("returned-->", repr(d))
    # print(`d[:]`)

    f.close()
