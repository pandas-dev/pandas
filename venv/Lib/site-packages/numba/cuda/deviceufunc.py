"""
Implements custom ufunc dispatch mechanism for non-CPU devices.
"""

from abc import ABCMeta, abstractmethod
from collections import OrderedDict
import operator
import warnings
from functools import reduce

import numpy as np

from numba.np.ufunc.ufuncbuilder import _BaseUFuncBuilder, parse_identity
from numba.core import types, sigutils
from numba.core.typing import signature
from numba.np.ufunc.sigparse import parse_signature


def _broadcast_axis(a, b):
    """
    Raises
    ------
    ValueError if broadcast fails
    """
    if a == b:
        return a
    elif a == 1:
        return b
    elif b == 1:
        return a
    else:
        raise ValueError("failed to broadcast {0} and {1}".format(a, b))


def _pairwise_broadcast(shape1, shape2):
    """
    Raises
    ------
    ValueError if broadcast fails
    """
    shape1, shape2 = map(tuple, [shape1, shape2])

    while len(shape1) < len(shape2):
        shape1 = (1,) + shape1

    while len(shape1) > len(shape2):
        shape2 = (1,) + shape2

    return tuple(_broadcast_axis(a, b) for a, b in zip(shape1, shape2))


def _multi_broadcast(*shapelist):
    """
    Raises
    ------
    ValueError if broadcast fails
    """
    assert shapelist

    result = shapelist[0]
    others = shapelist[1:]
    try:
        for i, each in enumerate(others, start=1):
            result = _pairwise_broadcast(result, each)
    except ValueError:
        raise ValueError("failed to broadcast argument #{0}".format(i))
    else:
        return result


class UFuncMechanism(object):
    """
    Prepare ufunc arguments for vectorize.
    """
    DEFAULT_STREAM = None
    SUPPORT_DEVICE_SLICING = False

    def __init__(self, typemap, args):
        """Never used directly by user. Invoke by UFuncMechanism.call().
        """
        self.typemap = typemap
        self.args = args
        nargs = len(self.args)
        self.argtypes = [None] * nargs
        self.scalarpos = []
        self.signature = None
        self.arrays = [None] * nargs

    def _fill_arrays(self):
        """
        Get all arguments in array form
        """
        for i, arg in enumerate(self.args):
            if self.is_device_array(arg):
                self.arrays[i] = self.as_device_array(arg)
            elif isinstance(arg, (int, float, complex, np.number)):
                # Is scalar
                self.scalarpos.append(i)
            else:
                self.arrays[i] = np.asarray(arg)

    def _fill_argtypes(self):
        """
        Get dtypes
        """
        for i, ary in enumerate(self.arrays):
            if ary is not None:
                dtype = getattr(ary, 'dtype')
                if dtype is None:
                    dtype = np.asarray(ary).dtype
                self.argtypes[i] = dtype

    def _resolve_signature(self):
        """Resolve signature.
        May have ambiguous case.
        """
        matches = []
        # Resolve scalar args exact match first
        if self.scalarpos:
            # Try resolve scalar arguments
            for formaltys in self.typemap:
                match_map = []
                for i, (formal, actual) in enumerate(zip(formaltys,
                                                         self.argtypes)):
                    if actual is None:
                        actual = np.asarray(self.args[i]).dtype

                    match_map.append(actual == formal)

                if all(match_map):
                    matches.append(formaltys)

        # No matching with exact match; try coercing the scalar arguments
        if not matches:
            matches = []
            for formaltys in self.typemap:
                all_matches = all(actual is None or formal == actual
                                  for formal, actual in
                                  zip(formaltys, self.argtypes))
                if all_matches:
                    matches.append(formaltys)

        if not matches:
            raise TypeError("No matching version.  GPU ufunc requires array "
                            "arguments to have the exact types.  This behaves "
                            "like regular ufunc with casting='no'.")

        if len(matches) > 1:
            raise TypeError("Failed to resolve ufunc due to ambiguous "
                            "signature. Too many untyped scalars. "
                            "Use numpy dtype object to type tag.")

        # Try scalar arguments
        self.argtypes = matches[0]

    def _get_actual_args(self):
        """Return the actual arguments
        Casts scalar arguments to np.array.
        """
        for i in self.scalarpos:
            self.arrays[i] = np.array([self.args[i]], dtype=self.argtypes[i])

        return self.arrays

    def _broadcast(self, arys):
        """Perform numpy ufunc broadcasting
        """
        shapelist = [a.shape for a in arys]
        shape = _multi_broadcast(*shapelist)

        for i, ary in enumerate(arys):
            if ary.shape == shape:
                pass

            else:
                if self.is_device_array(ary):
                    arys[i] = self.broadcast_device(ary, shape)

                else:
                    ax_differs = [ax for ax in range(len(shape))
                                  if ax >= ary.ndim
                                  or ary.shape[ax] != shape[ax]]

                    missingdim = len(shape) - len(ary.shape)
                    strides = [0] * missingdim + list(ary.strides)

                    for ax in ax_differs:
                        strides[ax] = 0

                    strided = np.lib.stride_tricks.as_strided(ary,
                                                              shape=shape,
                                                              strides=strides)

                    arys[i] = self.force_array_layout(strided)

        return arys

    def get_arguments(self):
        """Prepare and return the arguments for the ufunc.
        Does not call to_device().
        """
        self._fill_arrays()
        self._fill_argtypes()
        self._resolve_signature()
        arys = self._get_actual_args()
        return self._broadcast(arys)

    def get_function(self):
        """Returns (result_dtype, function)
        """
        return self.typemap[self.argtypes]

    def is_device_array(self, obj):
        """Is the `obj` a device array?
        Override in subclass
        """
        return False

    def as_device_array(self, obj):
        """Convert the `obj` to a device array
        Override in subclass

        Default implementation is an identity function
        """
        return obj

    def broadcast_device(self, ary, shape):
        """Handles ondevice broadcasting

        Override in subclass to add support.
        """
        raise NotImplementedError("broadcasting on device is not supported")

    def force_array_layout(self, ary):
        """Ensures array layout met device requirement.

        Override in sublcass
        """
        return ary

    @classmethod
    def call(cls, typemap, args, kws):
        """Perform the entire ufunc call mechanism.
        """
        # Handle keywords
        stream = kws.pop('stream', cls.DEFAULT_STREAM)
        out = kws.pop('out', None)

        if kws:
            warnings.warn("unrecognized keywords: %s" % ', '.join(kws))

        # Begin call resolution
        cr = cls(typemap, args)
        args = cr.get_arguments()
        resty, func = cr.get_function()

        outshape = args[0].shape

        # Adjust output value
        if out is not None and cr.is_device_array(out):
            out = cr.as_device_array(out)

        def attempt_ravel(a):
            if cr.SUPPORT_DEVICE_SLICING:
                raise NotImplementedError

            try:
                # Call the `.ravel()` method
                return a.ravel()
            except NotImplementedError:
                # If it is not a device array
                if not cr.is_device_array(a):
                    raise
                # For device array, retry ravel on the host by first
                # copying it back.
                else:
                    hostary = cr.to_host(a, stream).ravel()
                    return cr.to_device(hostary, stream)

        if args[0].ndim > 1:
            args = [attempt_ravel(a) for a in args]

        # Prepare argument on the device
        devarys = []
        any_device = False
        for a in args:
            if cr.is_device_array(a):
                devarys.append(a)
                any_device = True
            else:
                dev_a = cr.to_device(a, stream=stream)
                devarys.append(dev_a)

        # Launch
        shape = args[0].shape
        if out is None:
            # No output is provided
            devout = cr.allocate_device_array(shape, resty, stream=stream)

            devarys.extend([devout])
            cr.launch(func, shape[0], stream, devarys)

            if any_device:
                # If any of the arguments are on device,
                # Keep output on the device
                return devout.reshape(outshape)
            else:
                # Otherwise, transfer output back to host
                return devout.copy_to_host().reshape(outshape)

        elif cr.is_device_array(out):
            # If output is provided and it is a device array,
            # Return device array
            if out.ndim > 1:
                out = attempt_ravel(out)
            devout = out
            devarys.extend([devout])
            cr.launch(func, shape[0], stream, devarys)
            return devout.reshape(outshape)

        else:
            # If output is provided and it is a host array,
            # Return host array
            assert out.shape == shape
            assert out.dtype == resty
            devout = cr.allocate_device_array(shape, resty, stream=stream)
            devarys.extend([devout])
            cr.launch(func, shape[0], stream, devarys)
            return devout.copy_to_host(out, stream=stream).reshape(outshape)

    def to_device(self, hostary, stream):
        """Implement to device transfer
        Override in subclass
        """
        raise NotImplementedError

    def to_host(self, devary, stream):
        """Implement to host transfer
        Override in subclass
        """
        raise NotImplementedError

    def allocate_device_array(self, shape, dtype, stream):
        """Implements device allocation
        Override in subclass
        """
        raise NotImplementedError

    def launch(self, func, count, stream, args):
        """Implements device function invocation
        Override in subclass
        """
        raise NotImplementedError


def to_dtype(ty):
    if isinstance(ty, types.EnumMember):
        ty = ty.dtype
    return np.dtype(str(ty))


class DeviceVectorize(_BaseUFuncBuilder):
    def __init__(self, func, identity=None, cache=False, targetoptions=None):
        if targetoptions is None:
            targetoptions = {}
        if cache:
            raise TypeError("caching is not supported")
        for opt in targetoptions:
            if opt == 'nopython':
                warnings.warn("nopython kwarg for cuda target is redundant",
                              RuntimeWarning)
            else:
                fmt = "Unrecognized options. "
                fmt += "cuda vectorize target does not support option: '%s'"
                raise KeyError(fmt % opt)
        self.py_func = func
        self.identity = parse_identity(identity)
        # { arg_dtype: (return_dtype), cudakernel }
        self.kernelmap = OrderedDict()

    @property
    def pyfunc(self):
        return self.py_func

    def add(self, sig=None):
        # compile core as device function
        args, return_type = sigutils.normalize_signature(sig)
        devfnsig = signature(return_type, *args)

        funcname = self.pyfunc.__name__
        kernelsource = self._get_kernel_source(self._kernel_template,
                                               devfnsig, funcname)
        corefn, return_type = self._compile_core(devfnsig)
        glbl = self._get_globals(corefn)
        sig = signature(types.void, *([a[:] for a in args] + [return_type[:]]))
        exec(kernelsource, glbl)

        stager = glbl['__vectorized_%s' % funcname]
        kernel = self._compile_kernel(stager, sig)

        argdtypes = tuple(to_dtype(t) for t in devfnsig.args)
        resdtype = to_dtype(return_type)
        self.kernelmap[tuple(argdtypes)] = resdtype, kernel

    def build_ufunc(self):
        raise NotImplementedError

    def _get_kernel_source(self, template, sig, funcname):
        args = ['a%d' % i for i in range(len(sig.args))]
        fmts = dict(name=funcname,
                    args=', '.join(args),
                    argitems=', '.join('%s[__tid__]' % i for i in args))
        return template.format(**fmts)

    def _compile_core(self, sig):
        raise NotImplementedError

    def _get_globals(self, corefn):
        raise NotImplementedError

    def _compile_kernel(self, fnobj, sig):
        raise NotImplementedError


class DeviceGUFuncVectorize(_BaseUFuncBuilder):
    def __init__(
        self,
        func,
        sig,
        identity=None,
        cache=False,
        targetoptions=None,
        writable_args=(),
    ):
        if targetoptions is None:
            targetoptions = {}
        if cache:
            raise TypeError("caching is not supported")
        if writable_args:
            raise TypeError("writable_args are not supported")

        # Allow nopython flag to be set.
        if not targetoptions.pop('nopython', True):
            raise TypeError("nopython flag must be True")
        # Are there any more target options?
        if targetoptions:
            opts = ', '.join([repr(k) for k in targetoptions.keys()])
            fmt = "The following target options are not supported: {0}"
            raise TypeError(fmt.format(opts))

        self.py_func = func
        self.identity = parse_identity(identity)
        self.signature = sig
        self.inputsig, self.outputsig = parse_signature(self.signature)

        # Maps from a tuple of input_dtypes to (output_dtypes, kernel)
        self.kernelmap = OrderedDict()

    @property
    def pyfunc(self):
        return self.py_func

    def add(self, sig=None):
        indims = [len(x) for x in self.inputsig]
        outdims = [len(x) for x in self.outputsig]
        args, return_type = sigutils.normalize_signature(sig)

        # It is only valid to specify types.none as a return type, or to not
        # specify the return type (where the "Python None" is the return type)
        valid_return_type = return_type in (types.none, None)
        if not valid_return_type:
            raise TypeError('guvectorized functions cannot return values: '
                            f'signature {sig} specifies {return_type} return '
                            'type')

        funcname = self.py_func.__name__
        src = expand_gufunc_template(self._kernel_template, indims,
                                     outdims, funcname, args)

        glbls = self._get_globals(sig)

        exec(src, glbls)
        fnobj = glbls['__gufunc_{name}'.format(name=funcname)]

        outertys = list(_determine_gufunc_outer_types(args, indims + outdims))
        kernel = self._compile_kernel(fnobj, sig=tuple(outertys))

        nout = len(outdims)
        dtypes = [np.dtype(str(t.dtype)) for t in outertys]
        indtypes = tuple(dtypes[:-nout])
        outdtypes = tuple(dtypes[-nout:])

        self.kernelmap[indtypes] = outdtypes, kernel

    def _compile_kernel(self, fnobj, sig):
        raise NotImplementedError

    def _get_globals(self, sig):
        raise NotImplementedError


def _determine_gufunc_outer_types(argtys, dims):
    for at, nd in zip(argtys, dims):
        if isinstance(at, types.Array):
            yield at.copy(ndim=nd + 1)
        else:
            if nd > 0:
                raise ValueError("gufunc signature mismatch: ndim>0 for scalar")
            yield types.Array(dtype=at, ndim=1, layout='A')


def expand_gufunc_template(template, indims, outdims, funcname, argtypes):
    """Expand gufunc source template
    """
    argdims = indims + outdims
    argnames = ["arg{0}".format(i) for i in range(len(argdims))]
    checkedarg = "min({0})".format(', '.join(["{0}.shape[0]".format(a)
                                              for a in argnames]))
    inputs = [_gen_src_for_indexing(aref, adims, atype)
              for aref, adims, atype in zip(argnames, indims, argtypes)]
    outputs = [_gen_src_for_indexing(aref, adims, atype)
               for aref, adims, atype in zip(argnames[len(indims):], outdims,
                                             argtypes[len(indims):])]
    argitems = inputs + outputs
    src = template.format(name=funcname, args=', '.join(argnames),
                          checkedarg=checkedarg,
                          argitems=', '.join(argitems))
    return src


def _gen_src_for_indexing(aref, adims, atype):
    return "{aref}[{sliced}]".format(aref=aref,
                                     sliced=_gen_src_index(adims, atype))


def _gen_src_index(adims, atype):
    if adims > 0:
        return ','.join(['__tid__'] + [':'] * adims)
    elif isinstance(atype, types.Array) and atype.ndim - 1 == adims:
        # Special case for 0-nd in shape-signature but
        # 1d array in type signature.
        # Slice it so that the result has the same dimension.
        return '__tid__:(__tid__ + 1)'
    else:
        return '__tid__'


class GUFuncEngine(object):
    '''Determine how to broadcast and execute a gufunc
    base on input shape and signature
    '''

    @classmethod
    def from_signature(cls, signature):
        return cls(*parse_signature(signature))

    def __init__(self, inputsig, outputsig):
        # signatures
        self.sin = inputsig
        self.sout = outputsig
        # argument count
        self.nin = len(self.sin)
        self.nout = len(self.sout)

    def schedule(self, ishapes):
        if len(ishapes) != self.nin:
            raise TypeError('invalid number of input argument')

        # associate symbol values for input signature
        symbolmap = {}
        outer_shapes = []
        inner_shapes = []

        for argn, (shape, symbols) in enumerate(zip(ishapes, self.sin)):
            argn += 1  # start from 1 for human
            inner_ndim = len(symbols)
            if len(shape) < inner_ndim:
                fmt = "arg #%d: insufficient inner dimension"
                raise ValueError(fmt % (argn,))
            if inner_ndim:
                inner_shape = shape[-inner_ndim:]
                outer_shape = shape[:-inner_ndim]
            else:
                inner_shape = ()
                outer_shape = shape

            for axis, (dim, sym) in enumerate(zip(inner_shape, symbols)):
                axis += len(outer_shape)
                if sym in symbolmap:
                    if symbolmap[sym] != dim:
                        fmt = "arg #%d: shape[%d] mismatch argument"
                        raise ValueError(fmt % (argn, axis))
                symbolmap[sym] = dim

            outer_shapes.append(outer_shape)
            inner_shapes.append(inner_shape)

        # solve output shape
        oshapes = []
        for outsig in self.sout:
            oshape = []
            for sym in outsig:
                oshape.append(symbolmap[sym])
            oshapes.append(tuple(oshape))

        # find the biggest outershape as looping dimension
        sizes = [reduce(operator.mul, s, 1) for s in outer_shapes]
        largest_i = np.argmax(sizes)
        loopdims = outer_shapes[largest_i]

        pinned = [False] * self.nin  # same argument for each iteration
        for i, d in enumerate(outer_shapes):
            if d != loopdims:
                if d == (1,) or d == ():
                    pinned[i] = True
                else:
                    fmt = "arg #%d: outer dimension mismatch"
                    raise ValueError(fmt % (i + 1,))

        return GUFuncSchedule(self, inner_shapes, oshapes, loopdims, pinned)


class GUFuncSchedule(object):
    def __init__(self, parent, ishapes, oshapes, loopdims, pinned):
        self.parent = parent
        # core shapes
        self.ishapes = ishapes
        self.oshapes = oshapes
        # looping dimension
        self.loopdims = loopdims
        self.loopn = reduce(operator.mul, loopdims, 1)
        # flags
        self.pinned = pinned

        self.output_shapes = [loopdims + s for s in oshapes]

    def __str__(self):
        import pprint

        attrs = 'ishapes', 'oshapes', 'loopdims', 'loopn', 'pinned'
        values = [(k, getattr(self, k)) for k in attrs]
        return pprint.pformat(dict(values))


class GeneralizedUFunc(object):
    def __init__(self, kernelmap, engine):
        self.kernelmap = kernelmap
        self.engine = engine
        self.max_blocksize = 2 ** 30

    def __call__(self, *args, **kws):
        callsteps = self._call_steps(self.engine.nin, self.engine.nout,
                                     args, kws)
        indtypes, schedule, outdtypes, kernel = self._schedule(
            callsteps.inputs, callsteps.outputs)
        callsteps.adjust_input_types(indtypes)

        outputs = callsteps.prepare_outputs(schedule, outdtypes)
        inputs = callsteps.prepare_inputs()
        parameters = self._broadcast(schedule, inputs, outputs)

        callsteps.launch_kernel(kernel, schedule.loopn, parameters)

        return callsteps.post_process_outputs(outputs)

    def _schedule(self, inputs, outs):
        input_shapes = [a.shape for a in inputs]
        schedule = self.engine.schedule(input_shapes)

        # find kernel
        indtypes = tuple(i.dtype for i in inputs)
        try:
            outdtypes, kernel = self.kernelmap[indtypes]
        except KeyError:
            # No exact match, then use the first compatible.
            # This does not match the numpy dispatching exactly.
            # Later, we may just jit a new version for the missing signature.
            indtypes = self._search_matching_signature(indtypes)
            # Select kernel
            outdtypes, kernel = self.kernelmap[indtypes]

        # check output
        for sched_shape, out in zip(schedule.output_shapes, outs):
            if out is not None and sched_shape != out.shape:
                raise ValueError('output shape mismatch')

        return indtypes, schedule, outdtypes, kernel

    def _search_matching_signature(self, idtypes):
        """
        Given the input types in `idtypes`, return a compatible sequence of
        types that is defined in `kernelmap`.

        Note: Ordering is guaranteed by `kernelmap` being a OrderedDict
        """
        for sig in self.kernelmap.keys():
            if all(np.can_cast(actual, desired)
                   for actual, desired in zip(sig, idtypes)):
                return sig
        else:
            raise TypeError("no matching signature")

    def _broadcast(self, schedule, params, retvals):
        assert schedule.loopn > 0, "zero looping dimension"

        odim = 1 if not schedule.loopdims else schedule.loopn
        newparams = []
        for p, cs in zip(params, schedule.ishapes):
            if not cs and p.size == 1:
                # Broadcast scalar input
                devary = self._broadcast_scalar_input(p, odim)
                newparams.append(devary)
            else:
                # Broadcast vector input
                newparams.append(self._broadcast_array(p, odim, cs))

        newretvals = []
        for retval, oshape in zip(retvals, schedule.oshapes):
            newretvals.append(retval.reshape(odim, *oshape))
        return tuple(newparams) + tuple(newretvals)

    def _broadcast_array(self, ary, newdim, innerdim):
        newshape = (newdim,) + innerdim
        # No change in shape
        if ary.shape == newshape:
            return ary

        # Creating new dimension
        elif len(ary.shape) < len(newshape):
            assert newshape[-len(ary.shape):] == ary.shape, \
                "cannot add dim and reshape at the same time"
            return self._broadcast_add_axis(ary, newshape)

        # Collapsing dimension
        else:
            return ary.reshape(*newshape)

    def _broadcast_add_axis(self, ary, newshape):
        raise NotImplementedError("cannot add new axis")

    def _broadcast_scalar_input(self, ary, shape):
        raise NotImplementedError


class GUFuncCallSteps(metaclass=ABCMeta):
    """
    Implements memory management and kernel launch operations for GUFunc calls.

    One instance of this class is instantiated for each call, and the instance
    is specific to the arguments given to the GUFunc call.

    The base class implements the overall logic; subclasses provide
    target-specific implementations of individual functions.
    """

    # The base class uses these slots; subclasses may provide additional slots.
    __slots__ = [
        'outputs',
        'inputs',
        '_copy_result_to_host',
    ]

    @abstractmethod
    def launch_kernel(self, kernel, nelem, args):
        """Implement the kernel launch"""

    @abstractmethod
    def is_device_array(self, obj):
        """
        Return True if `obj` is a device array for this target, False
        otherwise.
        """

    @abstractmethod
    def as_device_array(self, obj):
        """
        Return `obj` as a device array on this target.

        May return `obj` directly if it is already on the target.
        """

    @abstractmethod
    def to_device(self, hostary):
        """
        Copy `hostary` to the device and return the device array.
        """

    @abstractmethod
    def allocate_device_array(self, shape, dtype):
        """
        Allocate a new uninitialized device array with the given shape and
        dtype.
        """

    def __init__(self, nin, nout, args, kwargs):
        outputs = kwargs.get('out')

        # Ensure the user has passed a correct number of arguments
        if outputs is None and len(args) not in (nin, (nin + nout)):
            def pos_argn(n):
                return f'{n} positional argument{"s" * (n != 1)}'

            msg = (f'This gufunc accepts {pos_argn(nin)} (when providing '
                   f'input only) or {pos_argn(nin + nout)} (when providing '
                   f'input and output). Got {pos_argn(len(args))}.')
            raise TypeError(msg)

        if outputs is not None and len(args) > nin:
            raise ValueError("cannot specify argument 'out' as both positional "
                             "and keyword")
        else:
            # If the user did not pass outputs either in the out kwarg or as
            # positional arguments, then we need to generate an initial list of
            # "placeholder" outputs using None as a sentry value
            outputs = [outputs] * nout

        # Ensure all output device arrays are Numba device arrays - for
        # example, any output passed in that supports the CUDA Array Interface
        # is converted to a Numba CUDA device array; others are left untouched.
        all_user_outputs_are_host = True
        self.outputs = []
        for output in outputs:
            if self.is_device_array(output):
                self.outputs.append(self.as_device_array(output))
                all_user_outputs_are_host = False
            else:
                self.outputs.append(output)

        all_host_arrays = not any([self.is_device_array(a) for a in args])

        # - If any of the arguments are device arrays, we leave the output on
        #   the device.
        self._copy_result_to_host = (all_host_arrays and
                                     all_user_outputs_are_host)

        # Normalize arguments - ensure they are either device- or host-side
        # arrays (as opposed to lists, tuples, etc).
        def normalize_arg(a):
            if self.is_device_array(a):
                convert = self.as_device_array
            else:
                convert = np.asarray

            return convert(a)

        normalized_args = [normalize_arg(a) for a in args]
        self.inputs = normalized_args[:nin]

        # Check if there are extra arguments for outputs.
        unused_inputs = normalized_args[nin:]
        if unused_inputs:
            self.outputs = unused_inputs

    def adjust_input_types(self, indtypes):
        """
        Attempt to cast the inputs to the required types if necessary
        and if they are not device arrays.

        Side effect: Only affects the elements of `inputs` that require
        a type cast.
        """
        for i, (ity, val) in enumerate(zip(indtypes, self.inputs)):
            if ity != val.dtype:
                if not hasattr(val, 'astype'):
                    msg = ("compatible signature is possible by casting but "
                           "{0} does not support .astype()").format(type(val))
                    raise TypeError(msg)
                # Cast types
                self.inputs[i] = val.astype(ity)

    def prepare_outputs(self, schedule, outdtypes):
        """
        Returns a list of output parameters that all reside on the target
        device.

        Outputs that were passed-in to the GUFunc are used if they reside on the
        device; other outputs are allocated as necessary.
        """
        outputs = []
        for shape, dtype, output in zip(schedule.output_shapes, outdtypes,
                                        self.outputs):
            if output is None or self._copy_result_to_host:
                output = self.allocate_device_array(shape, dtype)
            outputs.append(output)

        return outputs

    def prepare_inputs(self):
        """
        Returns a list of input parameters that all reside on the target device.
        """
        def ensure_device(parameter):
            if self.is_device_array(parameter):
                convert = self.as_device_array
            else:
                convert = self.to_device

            return convert(parameter)

        return [ensure_device(p) for p in self.inputs]

    def post_process_outputs(self, outputs):
        """
        Moves the given output(s) to the host if necessary.

        Returns a single value (e.g. an array) if there was one output, or a
        tuple of arrays if there were multiple. Although this feels a little
        jarring, it is consistent with the behavior of GUFuncs in general.
        """
        if self._copy_result_to_host:
            outputs = [self.to_host(output, self_output)
                       for output, self_output in zip(outputs, self.outputs)]
        elif self.outputs[0] is not None:
            outputs = self.outputs

        if len(outputs) == 1:
            return outputs[0]
        else:
            return tuple(outputs)
