from numba.extending import (models, register_model, type_callable,
                             unbox, NativeValue, make_attribute_wrapper, box,
                             lower_builtin)
from numba.core import types, cgutils
import warnings
from numba.core.errors import NumbaExperimentalFeatureWarning, NumbaValueError
from numpy.polynomial.polynomial import Polynomial
from contextlib import ExitStack
import numpy as np
from llvmlite import ir


@register_model(types.PolynomialType)
class PolynomialModel(models.StructModel):
    def __init__(self, dmm, fe_type):
        members = [
            ('coef', fe_type.coef),
            ('domain', fe_type.domain),
            ('window', fe_type.window)
            # Introduced in NumPy 1.24, maybe leave it out for now
            # ('symbol', types.string)
        ]
        super(PolynomialModel, self).__init__(dmm, fe_type, members)


@type_callable(Polynomial)
def type_polynomial(context):
    def typer(coef, domain=None, window=None):
        default_domain = types.Array(types.int64, 1, 'C')
        double_domain = types.Array(types.double, 1, 'C')
        default_window = types.Array(types.int64, 1, 'C')
        double_window = types.Array(types.double, 1, 'C')
        double_coef = types.Array(types.double, 1, 'C')

        warnings.warn("Polynomial class is experimental",
                      category=NumbaExperimentalFeatureWarning)

        if isinstance(coef, types.Array) and \
                all([a is None for a in (domain, window)]):
            if coef.ndim == 1:
                # If Polynomial(coef) is called, coef is cast to double dtype,
                # and domain and window are set to equal [-1, 1], i.e. have
                # integer dtype
                return types.PolynomialType(double_coef,
                                            default_domain,
                                            default_window,
                                            1)
            else:
                msg = 'Coefficient array is not 1-d'
                raise NumbaValueError(msg)
        elif all([isinstance(a, types.Array) for a in (coef, domain, window)]):
            if coef.ndim == 1:
                if all([a.ndim == 1 for a in (domain, window)]):
                    # If Polynomial(coef, domain, window) is called, then coef,
                    # domain and window are cast to double dtype
                    return types.PolynomialType(double_coef,
                                                double_domain,
                                                double_window,
                                                3)
            else:
                msg = 'Coefficient array is not 1-d'
                raise NumbaValueError(msg)
    return typer


make_attribute_wrapper(types.PolynomialType, 'coef', 'coef')
make_attribute_wrapper(types.PolynomialType, 'domain', 'domain')
make_attribute_wrapper(types.PolynomialType, 'window', 'window')
# Introduced in NumPy 1.24, maybe leave it out for now
# make_attribute_wrapper(types.PolynomialType, 'symbol', 'symbol')


@lower_builtin(Polynomial, types.Array)
def impl_polynomial1(context, builder, sig, args):

    def to_double(arr):
        return np.asarray(arr, dtype=np.double)

    def const_impl():
        return np.asarray([-1, 1])

    typ = sig.return_type
    polynomial = cgutils.create_struct_proxy(typ)(context, builder)
    sig_coef = sig.args[0].copy(dtype=types.double)(sig.args[0])
    coef_cast = context.compile_internal(builder, to_double, sig_coef, args)
    sig_domain = sig.args[0].copy(dtype=types.intp)()
    sig_window = sig.args[0].copy(dtype=types.intp)()
    domain_cast = context.compile_internal(builder, const_impl, sig_domain, ())
    window_cast = context.compile_internal(builder, const_impl, sig_window, ())
    polynomial.coef = coef_cast
    polynomial.domain = domain_cast
    polynomial.window = window_cast

    return polynomial._getvalue()


@lower_builtin(Polynomial, types.Array, types.Array, types.Array)
def impl_polynomial3(context, builder, sig, args):

    def to_double(coef):
        return np.asarray(coef, dtype=np.double)

    typ = sig.return_type
    polynomial = cgutils.create_struct_proxy(typ)(context, builder)

    coef_sig = sig.args[0].copy(dtype=types.double)(sig.args[0])
    domain_sig = sig.args[1].copy(dtype=types.double)(sig.args[1])
    window_sig = sig.args[2].copy(dtype=types.double)(sig.args[2])
    coef_cast = context.compile_internal(builder,
                                         to_double, coef_sig,
                                         (args[0],))
    domain_cast = context.compile_internal(builder,
                                           to_double, domain_sig,
                                           (args[1],))
    window_cast = context.compile_internal(builder,
                                           to_double, window_sig,
                                           (args[2],))

    domain_helper = context.make_helper(builder,
                                        domain_sig.return_type,
                                        value=domain_cast)
    window_helper = context.make_helper(builder,
                                        window_sig.return_type,
                                        value=window_cast)

    i64 = ir.IntType(64)
    two = i64(2)

    s1 = builder.extract_value(domain_helper.shape, 0)
    s2 = builder.extract_value(window_helper.shape, 0)
    pred1 = builder.icmp_signed('!=', s1, two)
    pred2 = builder.icmp_signed('!=', s2, two)

    with cgutils.if_unlikely(builder, pred1):
        context.call_conv.return_user_exc(
            builder, ValueError,
            ("Domain has wrong number of elements.",))

    with cgutils.if_unlikely(builder, pred2):
        context.call_conv.return_user_exc(
            builder, ValueError,
            ("Window has wrong number of elements.",))

    polynomial.coef = coef_cast
    polynomial.domain = domain_helper._getvalue()
    polynomial.window = window_helper._getvalue()

    return polynomial._getvalue()


@unbox(types.PolynomialType)
def unbox_polynomial(typ, obj, c):
    """
    Convert a Polynomial object to a native polynomial structure.
    """
    is_error_ptr = cgutils.alloca_once_value(c.builder, cgutils.false_bit)
    polynomial = cgutils.create_struct_proxy(typ)(c.context, c.builder)
    with ExitStack() as stack:
        natives = []
        for name in ("coef", "domain", "window"):
            attr = c.pyapi.object_getattr_string(obj, name)
            with cgutils.early_exit_if_null(c.builder, stack, attr):
                c.builder.store(cgutils.true_bit, is_error_ptr)
            t = getattr(typ, name)
            native = c.unbox(t, attr)
            c.pyapi.decref(attr)
            with cgutils.early_exit_if(c.builder, stack, native.is_error):
                c.builder.store(cgutils.true_bit, is_error_ptr)
            natives.append(native)

        polynomial.coef = natives[0]
        polynomial.domain = natives[1]
        polynomial.window = natives[2]

    return NativeValue(polynomial._getvalue(),
                       is_error=c.builder.load(is_error_ptr))


@box(types.PolynomialType)
def box_polynomial(typ, val, c):
    """
    Convert a native polynomial structure to a Polynomial object.
    """
    ret_ptr = cgutils.alloca_once(c.builder, c.pyapi.pyobj)
    fail_obj = c.pyapi.get_null_object()

    with ExitStack() as stack:
        polynomial = cgutils.create_struct_proxy(typ)(c.context, c.builder,
                                                      value=val)
        coef_obj = c.box(typ.coef, polynomial.coef)
        with cgutils.early_exit_if_null(c.builder, stack, coef_obj):
            c.builder.store(fail_obj, ret_ptr)

        domain_obj = c.box(typ.domain, polynomial.domain)
        with cgutils.early_exit_if_null(c.builder, stack, domain_obj):
            c.builder.store(fail_obj, ret_ptr)

        window_obj = c.box(typ.window, polynomial.window)
        with cgutils.early_exit_if_null(c.builder, stack, window_obj):
            c.builder.store(fail_obj, ret_ptr)

        class_obj = c.pyapi.unserialize(c.pyapi.serialize_object(Polynomial))
        with cgutils.early_exit_if_null(c.builder, stack, class_obj):
            c.pyapi.decref(coef_obj)
            c.pyapi.decref(domain_obj)
            c.pyapi.decref(window_obj)
            c.builder.store(fail_obj, ret_ptr)

        if typ.n_args == 1:
            res1 = c.pyapi.call_function_objargs(class_obj, (coef_obj,))
            c.builder.store(res1, ret_ptr)
        else:
            res3 = c.pyapi.call_function_objargs(class_obj, (coef_obj,
                                                             domain_obj,
                                                             window_obj))
            c.builder.store(res3, ret_ptr)

        c.pyapi.decref(coef_obj)
        c.pyapi.decref(domain_obj)
        c.pyapi.decref(window_obj)
        c.pyapi.decref(class_obj)

    return c.builder.load(ret_ptr)
