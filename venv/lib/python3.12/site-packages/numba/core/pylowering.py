"""
Lowering implementation for object mode.
"""


import builtins
import operator
import inspect
from functools import cached_property

import llvmlite.ir

from numba.core import types, utils, ir, generators, cgutils
from numba.core.errors import (ForbiddenConstruct, LoweringError,
                               NumbaNotImplementedError)
from numba.core.lowering import BaseLower


# Issue #475: locals() is unsupported as calling it naively would give
# out wrong results.
_unsupported_builtins = set([locals])


class _Undefined:
    """
    A sentinel value for undefined variable created by Expr.undef.
    """
    def __repr__(self):
        return "<undefined>"


_UNDEFINED = _Undefined()


# Map operators to methods on the PythonAPI class
PYTHON_BINOPMAP = {
    operator.add: ("number_add", False),
    operator.sub: ("number_subtract", False),
    operator.mul: ("number_multiply", False),
    operator.truediv: ("number_truedivide", False),
    operator.floordiv: ("number_floordivide", False),
    operator.mod: ("number_remainder", False),
    operator.pow: ("number_power", False),
    operator.lshift: ("number_lshift", False),
    operator.rshift: ("number_rshift", False),
    operator.and_: ("number_and", False),
    operator.or_: ("number_or", False),
    operator.xor: ("number_xor", False),
    # inplace operators
    operator.iadd: ("number_add", True),
    operator.isub: ("number_subtract", True),
    operator.imul: ("number_multiply", True),
    operator.itruediv: ("number_truedivide", True),
    operator.ifloordiv: ("number_floordivide", True),
    operator.imod: ("number_remainder", True),
    operator.ipow: ("number_power", True),
    operator.ilshift: ("number_lshift", True),
    operator.irshift: ("number_rshift", True),
    operator.iand: ("number_and", True),
    operator.ior: ("number_or", True),
    operator.ixor: ("number_xor", True),
}

PYTHON_BINOPMAP[operator.matmul] = ("number_matrix_multiply", False)
PYTHON_BINOPMAP[operator.imatmul] = ("number_matrix_multiply", True)

PYTHON_COMPAREOPMAP = {
    operator.eq: '==',
    operator.ne: '!=',
    operator.lt: '<',
    operator.le: '<=',
    operator.gt: '>',
    operator.ge: '>=',
    operator.is_: 'is',
    operator.is_not: 'is not',
    operator.contains: 'in'
}

class PyLower(BaseLower):

    GeneratorLower = generators.PyGeneratorLower

    def init(self):
        # Strings to be frozen into the Environment object
        self._frozen_strings = set()

        self._live_vars = set()

    def pre_lower(self):
        super(PyLower, self).pre_lower()
        self.init_pyapi()

    def post_lower(self):
        pass

    def pre_block(self, block):
        self.init_vars(block)

    def lower_inst(self, inst):
        if isinstance(inst, ir.Assign):
            value = self.lower_assign(inst)
            self.storevar(value, inst.target.name)

        elif isinstance(inst, ir.SetItem):
            target = self.loadvar(inst.target.name)
            index = self.loadvar(inst.index.name)
            value = self.loadvar(inst.value.name)
            ok = self.pyapi.object_setitem(target, index, value)
            self.check_int_status(ok)

        elif isinstance(inst, ir.DelItem):
            target = self.loadvar(inst.target.name)
            index = self.loadvar(inst.index.name)
            ok = self.pyapi.object_delitem(target, index)
            self.check_int_status(ok)

        elif isinstance(inst, ir.SetAttr):
            target = self.loadvar(inst.target.name)
            value = self.loadvar(inst.value.name)
            ok = self.pyapi.object_setattr(target,
                                           self._freeze_string(inst.attr),
                                           value)
            self.check_int_status(ok)

        elif isinstance(inst, ir.DelAttr):
            target = self.loadvar(inst.target.name)
            ok = self.pyapi.object_delattr(target,
                                           self._freeze_string(inst.attr))
            self.check_int_status(ok)

        elif isinstance(inst, ir.StoreMap):
            dct = self.loadvar(inst.dct.name)
            key = self.loadvar(inst.key.name)
            value = self.loadvar(inst.value.name)
            ok = self.pyapi.dict_setitem(dct, key, value)
            self.check_int_status(ok)

        elif isinstance(inst, ir.Return):
            retval = self.loadvar(inst.value.name)
            if self.generator_info:
                # StopIteration
                # We own a reference to the "return value", but we
                # don't return it.
                self.pyapi.decref(retval)
                self.genlower.return_from_generator(self)
                return
            # No need to incref() as the reference is already owned.
            self.call_conv.return_value(self.builder, retval)

        elif isinstance(inst, ir.Branch):
            cond = self.loadvar(inst.cond.name)
            if cond.type == llvmlite.ir.IntType(1):
                istrue = cond
            else:
                istrue = self.pyapi.object_istrue(cond)
            zero = llvmlite.ir.Constant(istrue.type, None)
            pred = self.builder.icmp_unsigned('!=', istrue, zero)
            tr = self.blkmap[inst.truebr]
            fl = self.blkmap[inst.falsebr]
            self.builder.cbranch(pred, tr, fl)

        elif isinstance(inst, ir.Jump):
            target = self.blkmap[inst.target]
            self.builder.branch(target)

        elif isinstance(inst, ir.Del):
            self.delvar(inst.value)

        elif isinstance(inst, ir.PopBlock):
            pass # this is just a marker

        elif isinstance(inst, ir.Raise):
            if inst.exception is not None:
                exc = self.loadvar(inst.exception.name)
                # A reference will be stolen by raise_object() and another
                # by return_exception_raised().
                self.incref(exc)
            else:
                exc = None
            self.pyapi.raise_object(exc)
            self.return_exception_raised()

        else:
            msg = f"{type(inst)}, {inst}"
            raise NumbaNotImplementedError(msg)

    @cached_property
    def _omitted_typobj(self):
        """Return a `OmittedArg` type instance as a LLVM value suitable for
        testing at runtime.
        """
        from numba.core.dispatcher import OmittedArg
        return self.pyapi.unserialize(
            self.pyapi.serialize_object(OmittedArg))

    def lower_assign(self, inst):
        """
        The returned object must have a new reference
        """
        value = inst.value
        if isinstance(value, (ir.Const, ir.FreeVar)):
            return self.lower_const(value.value)
        elif isinstance(value, ir.Var):
            val = self.loadvar(value.name)
            self.incref(val)
            return val
        elif isinstance(value, ir.Expr):
            return self.lower_expr(value)
        elif isinstance(value, ir.Global):
            return self.lower_global(value.name, value.value)
        elif isinstance(value, ir.Yield):
            return self.lower_yield(value)
        elif isinstance(value, ir.Arg):
            param = self.func_ir.func_id.pysig.parameters.get(value.name)

            obj = self.fnargs[value.index]
            slot = cgutils.alloca_once_value(self.builder, obj)
            # Don't check for OmittedArg unless the argument has a default
            if param is not None and param.default is inspect.Parameter.empty:
                self.incref(obj)
                self.builder.store(obj, slot)
            else:
                # When an argument is omitted, the dispatcher hands it as
                # _OmittedArg(<default value>)
                typobj = self.pyapi.get_type(obj)
                is_omitted = self.builder.icmp_unsigned('==', typobj,
                                                        self._omitted_typobj)
                with self.builder.if_else(is_omitted, likely=False) as (omitted, present):
                    with present:
                        self.incref(obj)
                        self.builder.store(obj, slot)
                    with omitted:
                        # The argument is omitted => get the default value
                        obj = self.pyapi.object_getattr_string(obj, 'value')
                        self.builder.store(obj, slot)

            return self.builder.load(slot)
        else:
            raise NotImplementedError(type(value), value)

    def lower_yield(self, inst):
        yp = self.generator_info.yield_points[inst.index]
        assert yp.inst is inst
        self.genlower.init_generator_state(self)

        # Save live vars in state
        # We also need to save live vars that are del'ed afterwards.
        y = generators.LowerYield(self, yp, yp.live_vars | yp.weak_live_vars)
        y.lower_yield_suspend()
        # Yield to caller
        val = self.loadvar(inst.value.name)
        # Let caller own the reference
        self.pyapi.incref(val)
        self.call_conv.return_value(self.builder, val)

        # Resumption point
        y.lower_yield_resume()
        # None is returned by the yield expression
        return self.pyapi.make_none()

    def lower_binop(self, expr, op, inplace=False):
        lhs = self.loadvar(expr.lhs.name)
        rhs = self.loadvar(expr.rhs.name)
        assert not isinstance(op, str)
        if op in PYTHON_BINOPMAP:
            fname, inplace = PYTHON_BINOPMAP[op]
            fn = getattr(self.pyapi, fname)
            res = fn(lhs, rhs, inplace=inplace)
        else:
            # Assumed to be rich comparison
            fn = PYTHON_COMPAREOPMAP.get(expr.fn, expr.fn)
            if fn == 'in':      # 'in' and operator.contains have args reversed
                lhs, rhs = rhs, lhs
            res = self.pyapi.object_richcompare(lhs, rhs, fn)
        self.check_error(res)
        return res

    def lower_expr(self, expr):
        if expr.op == 'binop':
            return self.lower_binop(expr, expr.fn, inplace=False)
        elif expr.op == 'inplace_binop':
            return self.lower_binop(expr, expr.fn, inplace=True)
        elif expr.op == 'unary':
            value = self.loadvar(expr.value.name)
            if expr.fn == operator.neg:
                res = self.pyapi.number_negative(value)
            elif expr.fn == operator.pos:
                res = self.pyapi.number_positive(value)
            elif expr.fn == operator.not_:
                res = self.pyapi.object_not(value)
                self.check_int_status(res)
                res = self.pyapi.bool_from_bool(res)
            elif expr.fn == operator.invert:
                res = self.pyapi.number_invert(value)
            else:
                raise NotImplementedError(expr)
            self.check_error(res)
            return res
        elif expr.op == 'call':
            argvals = [self.loadvar(a.name) for a in expr.args]
            fn = self.loadvar(expr.func.name)
            args = self.pyapi.tuple_pack(argvals)
            if expr.vararg:
                # Expand *args
                varargs = self.pyapi.sequence_tuple(
                                self.loadvar(expr.vararg.name))
                new_args = self.pyapi.sequence_concat(args, varargs)
                self.decref(varargs)
                self.decref(args)
                args = new_args
            if not expr.kws:
                # No named arguments
                ret = self.pyapi.call(fn, args, None)
            else:
                # Named arguments
                keyvalues = [(k, self.loadvar(v.name)) for k, v in expr.kws]
                kws = self.pyapi.dict_pack(keyvalues)
                ret = self.pyapi.call(fn, args, kws)
                self.decref(kws)
            self.decref(args)
            self.check_error(ret)
            return ret
        elif expr.op == 'getattr':
            obj = self.loadvar(expr.value.name)
            res = self.pyapi.object_getattr(obj, self._freeze_string(expr.attr))
            self.check_error(res)
            return res
        elif expr.op == 'build_tuple':
            items = [self.loadvar(it.name) for it in expr.items]
            res = self.pyapi.tuple_pack(items)
            self.check_error(res)
            return res
        elif expr.op == 'build_list':
            items = [self.loadvar(it.name) for it in expr.items]
            res = self.pyapi.list_pack(items)
            self.check_error(res)
            return res
        elif expr.op == 'build_map':
            res = self.pyapi.dict_new(expr.size)
            self.check_error(res)
            for k, v in expr.items:
                key = self.loadvar(k.name)
                value = self.loadvar(v.name)
                ok = self.pyapi.dict_setitem(res, key, value)
                self.check_int_status(ok)
            return res
        elif expr.op == 'build_set':
            items = [self.loadvar(it.name) for it in expr.items]
            res = self.pyapi.set_new()
            self.check_error(res)
            for it in items:
                ok = self.pyapi.set_add(res, it)
                self.check_int_status(ok)
            return res
        elif expr.op == 'getiter':
            obj = self.loadvar(expr.value.name)
            res = self.pyapi.object_getiter(obj)
            self.check_error(res)
            return res
        elif expr.op == 'iternext':
            iterobj = self.loadvar(expr.value.name)
            item = self.pyapi.iter_next(iterobj)
            is_valid = cgutils.is_not_null(self.builder, item)
            pair = self.pyapi.tuple_new(2)
            with self.builder.if_else(is_valid) as (then, otherwise):
                with then:
                    self.pyapi.tuple_setitem(pair, 0, item)
                with otherwise:
                    self.check_occurred()
                    # Make the tuple valid by inserting None as dummy
                    # iteration "result" (it will be ignored).
                    self.pyapi.tuple_setitem(pair, 0, self.pyapi.make_none())
            self.pyapi.tuple_setitem(pair, 1, self.pyapi.bool_from_bool(is_valid))
            return pair
        elif expr.op == 'pair_first':
            pair = self.loadvar(expr.value.name)
            first = self.pyapi.tuple_getitem(pair, 0)
            self.incref(first)
            return first
        elif expr.op == 'pair_second':
            pair = self.loadvar(expr.value.name)
            second = self.pyapi.tuple_getitem(pair, 1)
            self.incref(second)
            return second
        elif expr.op == 'exhaust_iter':
            iterobj = self.loadvar(expr.value.name)
            tup = self.pyapi.sequence_tuple(iterobj)
            self.check_error(tup)
            # Check tuple size is as expected
            tup_size = self.pyapi.tuple_size(tup)
            expected_size = self.context.get_constant(types.intp, expr.count)
            has_wrong_size = self.builder.icmp_unsigned('!=',
                                               tup_size, expected_size)
            with cgutils.if_unlikely(self.builder, has_wrong_size):
                self.return_exception(ValueError)
            return tup
        elif expr.op == 'getitem':
            value = self.loadvar(expr.value.name)
            index = self.loadvar(expr.index.name)
            res = self.pyapi.object_getitem(value, index)
            self.check_error(res)
            return res
        elif expr.op == 'static_getitem':
            value = self.loadvar(expr.value.name)
            index = self.context.get_constant(types.intp, expr.index)
            indexobj = self.pyapi.long_from_ssize_t(index)
            self.check_error(indexobj)
            res = self.pyapi.object_getitem(value, indexobj)
            self.decref(indexobj)
            self.check_error(res)
            return res
        elif expr.op == 'getslice':
            target = self.loadvar(expr.target.name)
            start = self.loadvar(expr.start.name)
            stop = self.loadvar(expr.stop.name)

            slicefn = self.get_builtin_obj("slice")
            sliceobj = self.pyapi.call_function_objargs(slicefn, (start, stop))
            self.decref(slicefn)
            self.check_error(sliceobj)

            res = self.pyapi.object_getitem(target, sliceobj)
            self.check_error(res)

            return res

        elif expr.op == 'cast':
            val = self.loadvar(expr.value.name)
            self.incref(val)
            return val
        elif expr.op == 'phi':
            raise LoweringError("PHI not stripped")

        elif expr.op == 'null':
            # Make null value
            return cgutils.get_null_value(self.pyapi.pyobj)

        elif expr.op == 'undef':
            # Use a sentinel value for undefined variable
            return self.lower_const(_UNDEFINED)

        else:
            raise NotImplementedError(expr)

    def lower_const(self, const):
        # All constants are frozen inside the environment
        index = self.env_manager.add_const(const)
        ret = self.env_manager.read_const(index)
        self.check_error(ret)
        self.incref(ret)
        return ret

    def lower_global(self, name, value):
        """
        1) Check global scope dictionary.
        2) Check __builtins__.
            2a) is it a dictionary (for non __main__ module)
            2b) is it a module (for __main__ module)
        """
        moddict = self.get_module_dict()
        obj = self.pyapi.dict_getitem(moddict, self._freeze_string(name))
        self.incref(obj)  # obj is borrowed

        try:
            if value in _unsupported_builtins:
                raise ForbiddenConstruct("builtins %s() is not supported"
                                         % name, loc=self.loc)
        except TypeError:
            # `value` is unhashable, ignore
            pass

        if hasattr(builtins, name):
            obj_is_null = self.is_null(obj)
            bbelse = self.builder.basic_block

            with self.builder.if_then(obj_is_null):
                mod = self.pyapi.dict_getitem(moddict,
                                          self._freeze_string("__builtins__"))
                builtin = self.builtin_lookup(mod, name)
                bbif = self.builder.basic_block

            retval = self.builder.phi(self.pyapi.pyobj)
            retval.add_incoming(obj, bbelse)
            retval.add_incoming(builtin, bbif)

        else:
            retval = obj
            with cgutils.if_unlikely(self.builder, self.is_null(retval)):
                self.pyapi.raise_missing_global_error(name)
                self.return_exception_raised()

        return retval

    # -------------------------------------------------------------------------

    def get_module_dict(self):
        return self.env_body.globals

    def get_builtin_obj(self, name):
        # XXX The builtins dict could be bound into the environment
        moddict = self.get_module_dict()
        mod = self.pyapi.dict_getitem(moddict,
                                      self._freeze_string("__builtins__"))
        return self.builtin_lookup(mod, name)

    def builtin_lookup(self, mod, name):
        """
        Args
        ----
        mod:
            The __builtins__ dictionary or module, as looked up in
            a module's globals.
        name: str
            The object to lookup
        """
        fromdict = self.pyapi.dict_getitem(mod, self._freeze_string(name))
        self.incref(fromdict)       # fromdict is borrowed
        bbifdict = self.builder.basic_block

        with cgutils.if_unlikely(self.builder, self.is_null(fromdict)):
            # This happen if we are using the __main__ module
            frommod = self.pyapi.object_getattr(mod, self._freeze_string(name))

            with cgutils.if_unlikely(self.builder, self.is_null(frommod)):
                self.pyapi.raise_missing_global_error(name)
                self.return_exception_raised()

            bbifmod = self.builder.basic_block

        builtin = self.builder.phi(self.pyapi.pyobj)
        builtin.add_incoming(fromdict, bbifdict)
        builtin.add_incoming(frommod, bbifmod)

        return builtin

    def check_occurred(self):
        """
        Return if an exception occurred.
        """
        err_occurred = cgutils.is_not_null(self.builder,
                                           self.pyapi.err_occurred())

        with cgutils.if_unlikely(self.builder, err_occurred):
            self.return_exception_raised()

    def check_error(self, obj):
        """
        Return if *obj* is NULL.
        """
        with cgutils.if_unlikely(self.builder, self.is_null(obj)):
            self.return_exception_raised()

        return obj

    def check_int_status(self, num, ok_value=0):
        """
        Raise an exception if *num* is smaller than *ok_value*.
        """
        ok = llvmlite.ir.Constant(num.type, ok_value)
        pred = self.builder.icmp_signed('<', num, ok)
        with cgutils.if_unlikely(self.builder, pred):
            self.return_exception_raised()

    def is_null(self, obj):
        return cgutils.is_null(self.builder, obj)

    def return_exception_raised(self):
        """
        Return with the currently raised exception.
        """
        self.cleanup_vars()
        self.call_conv.return_exc(self.builder)

    def init_vars(self, block):
        """
        Initialize live variables for *block*.
        """
        self._live_vars = set(self.func_ir.get_block_entry_vars(block))

    def _getvar(self, name, ltype=None):
        if name not in self.varmap:
            self.varmap[name] = self.alloca(name, ltype=ltype)
        return self.varmap[name]

    def loadvar(self, name):
        """
        Load the llvm value of the variable named *name*.
        """
        # If this raises then the live variables analysis is wrong
        assert name in self._live_vars, name
        ptr = self.varmap[name]
        val = self.builder.load(ptr)
        with cgutils.if_unlikely(self.builder, self.is_null(val)):
            self.pyapi.raise_missing_name_error(name)
            self.return_exception_raised()
        return val

    def delvar(self, name):
        """
        Delete the variable slot with the given name. This will decref
        the corresponding Python object.
        """
        # If this raises then the live variables analysis is wrong
        self._live_vars.remove(name)
        ptr = self._getvar(name)  # initializes `name` if not already
        self.decref(self.builder.load(ptr))
        # This is a safety guard against double decref's, but really
        # the IR should be correct and have only one Del per variable
        # and code path.
        self.builder.store(cgutils.get_null_value(ptr.type.pointee), ptr)

    def storevar(self, value, name, clobber=False):
        """
        Stores a llvm value and allocate stack slot if necessary.
        The llvm value can be of arbitrary type.
        """
        is_redefine = name in self._live_vars and not clobber
        ptr = self._getvar(name, ltype=value.type)
        if is_redefine:
            old = self.builder.load(ptr)
        else:
            self._live_vars.add(name)
        assert value.type == ptr.type.pointee, (str(value.type),
                                                str(ptr.type.pointee))
        self.builder.store(value, ptr)
        # Safe to call decref even on non python object
        if is_redefine:
            self.decref(old)

    def cleanup_vars(self):
        """
        Cleanup live variables.
        """
        for name in self._live_vars:
            ptr = self._getvar(name)
            self.decref(self.builder.load(ptr))

    def alloca(self, name, ltype=None):
        """
        Allocate a stack slot and initialize it to NULL.
        The default is to allocate a pyobject pointer.
        Use ``ltype`` to override.
        """
        if ltype is None:
            ltype = self.context.get_value_type(types.pyobject)
        with self.builder.goto_block(self.entry_block):
            ptr = self.builder.alloca(ltype, name=name)
            self.builder.store(cgutils.get_null_value(ltype), ptr)
        return ptr

    def _alloca_var(self, name, fetype):
        # This is here for API compatibility with lowering.py::Lower.
        # NOTE: fetype is unused
        return self.alloca(name)

    def incref(self, value):
        self.pyapi.incref(value)

    def decref(self, value):
        """
        This is allow to be called on non pyobject pointer, in which case
        no code is inserted.
        """
        lpyobj = self.context.get_value_type(types.pyobject)
        if value.type == lpyobj:
            self.pyapi.decref(value)

    def _freeze_string(self, string):
        """
        Freeze a Python string object into the code.
        """
        return self.lower_const(string)
