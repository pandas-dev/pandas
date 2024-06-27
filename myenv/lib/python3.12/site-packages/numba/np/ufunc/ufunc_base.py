from numba.np import numpy_support
from numba.core import types


class UfuncLowererBase:
    '''Callable class responsible for lowering calls to a specific gufunc.
    '''
    def __init__(self, ufunc, make_kernel_fn, make_ufunc_kernel_fn):
        self.ufunc = ufunc
        self.make_ufunc_kernel_fn = make_ufunc_kernel_fn
        self.kernel = make_kernel_fn(ufunc)
        self.libs = []

    def __call__(self, context, builder, sig, args):
        return self.make_ufunc_kernel_fn(context, builder, sig, args,
                                         self.ufunc, self.kernel)


class UfuncBase:

    @property
    def nin(self):
        return self.ufunc.nin

    @property
    def nout(self):
        return self.ufunc.nout

    @property
    def nargs(self):
        return self.ufunc.nargs

    @property
    def ntypes(self):
        return self.ufunc.ntypes

    @property
    def types(self):
        return self.ufunc.types

    @property
    def identity(self):
        return self.ufunc.identity

    @property
    def signature(self):
        return self.ufunc.signature

    @property
    def accumulate(self):
        return self.ufunc.accumulate

    @property
    def at(self):
        return self.ufunc.at

    @property
    def outer(self):
        return self.ufunc.outer

    @property
    def reduce(self):
        return self.ufunc.reduce

    @property
    def reduceat(self):
        return self.ufunc.reduceat

    def disable_compile(self):
        """
        Disable the compilation of new signatures at call time.
        """
        # If disabling compilation then there must be at least one signature
        assert len(self._dispatcher.overloads) > 0
        self._frozen = True

    def _install_cg(self, targetctx=None):
        """
        Install an implementation function for a GUFunc/DUFunc object in the
        given target context.  If no target context is given, then
        _install_cg() installs into the target context of the
        dispatcher object (should be same default context used by
        jit() and njit()).
        """
        if targetctx is None:
            targetctx = self._dispatcher.targetdescr.target_context
        _any = types.Any
        _arr = types.Array
        # Either all outputs are explicit or none of them are
        sig0 = (_any,) * self.ufunc.nin + (_arr,) * self.ufunc.nout
        sig1 = (_any,) * self.ufunc.nin
        targetctx.insert_func_defn(
            [(self._lower_me, self, sig) for sig in (sig0, sig1)])

    def find_ewise_function(self, ewise_types):
        """
        Given a tuple of element-wise argument types, find a matching
        signature in the dispatcher.

        Return a 2-tuple containing the matching signature, and
        compilation result.  Will return two None's if no matching
        signature was found.
        """
        if self._frozen:
            # If we cannot compile, coerce to the best matching loop
            loop = numpy_support.ufunc_find_matching_loop(self, ewise_types)
            if loop is None:
                return None, None
            ewise_types = tuple(loop.inputs + loop.outputs)[:len(ewise_types)]
        for sig, cres in self._dispatcher.overloads.items():
            if self.match_signature(ewise_types, sig):
                return sig, cres
        return None, None
