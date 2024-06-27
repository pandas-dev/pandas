from types import ModuleType

import weakref

from numba.core.errors import ConstantInferenceError, NumbaError
from numba.core import ir


class ConstantInference(object):
    """
    A constant inference engine for a given interpreter.
    Inference inspects the IR to try and compute a compile-time constant for
    a variable.

    This shouldn't be used directly, instead call Interpreter.infer_constant().
    """

    def __init__(self, func_ir):
        # Avoid cyclic references as some user-visible objects may be
        # held alive in the cache
        self._func_ir = weakref.proxy(func_ir)
        self._cache = {}

    def infer_constant(self, name, loc=None):
        """
        Infer a constant value for the given variable *name*.
        If no value can be inferred, numba.errors.ConstantInferenceError
        is raised.
        """
        if name not in self._cache:
            try:
                self._cache[name] = (True, self._do_infer(name))
            except ConstantInferenceError as exc:
                # Store the exception args only, to avoid keeping
                # a whole traceback alive.
                self._cache[name] = (False, (exc.__class__, exc.args))
        success, val = self._cache[name]
        if success:
            return val
        else:
            exc, args = val
            if issubclass(exc, NumbaError):
                raise exc(*args, loc=loc)
            else:
                raise exc(*args)

    def _fail(self, val):
        # The location here is set to None because `val` is the ir.Var name
        # and not the actual offending use of the var. When this is raised it is
        # caught in the flow control of `infer_constant` and the class and args
        # (the message) are captured and then raised again but with the location
        # set to the expression that caused the constant inference error.
        raise ConstantInferenceError(
            "Constant inference not possible for: %s" % (val,), loc=None)

    def _do_infer(self, name):
        if not isinstance(name, str):
            raise TypeError("infer_constant() called with non-str %r"
                            % (name,))
        try:
            defn = self._func_ir.get_definition(name)
        except KeyError:
            raise ConstantInferenceError(
                "no single definition for %r" % (name,))
        try:
            const = defn.infer_constant()
        except ConstantInferenceError:
            if isinstance(defn, ir.Expr):
                return self._infer_expr(defn)
            self._fail(defn)
        return const

    def _infer_expr(self, expr):
        # Infer an expression: handle supported cases
        if expr.op == 'call':
            func = self.infer_constant(expr.func.name, loc=expr.loc)
            return self._infer_call(func, expr)
        elif expr.op == 'getattr':
            value = self.infer_constant(expr.value.name, loc=expr.loc)
            return self._infer_getattr(value, expr)
        elif expr.op == 'build_list':
            return [self.infer_constant(i.name, loc=expr.loc) for i in
                    expr.items]
        elif expr.op == 'build_tuple':
            return tuple(self.infer_constant(i.name, loc=expr.loc) for i in
                         expr.items)
        self._fail(expr)

    def _infer_call(self, func, expr):
        if expr.kws or expr.vararg:
            self._fail(expr)
        # Check supported callables
        _slice = func in (slice,)
        _exc = isinstance(func, type) and issubclass(func, BaseException)
        if _slice or _exc:
            args = [self.infer_constant(a.name, loc=expr.loc) for a in
                    expr.args]
            if _slice:
                return func(*args)
            elif _exc:
                # If the exception class is user defined it may implement a ctor
                # that does not pass the args to the super. Therefore return the
                # raw class and the args so this can be instantiated at the call
                # site in the way the user source expects it to be.
                return func, args
            else:
                assert 0, 'Unreachable'

        self._fail(expr)

    def _infer_getattr(self, value, expr):
        if isinstance(value, (ModuleType, type)):
            # Allow looking up a constant on a class or module
            try:
                return getattr(value, expr.attr)
            except AttributeError:
                pass
        self._fail(expr)
