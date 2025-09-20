from .. import types, utils, errors
import operator
from .templates import (AttributeTemplate, ConcreteTemplate, AbstractTemplate,
                        infer_global, infer, infer_getattr,
                        signature, bound_function, make_callable_template)
from .builtins import normalize_1d_index


@infer_global(operator.contains)
class InContainer(AbstractTemplate):
    key = operator.contains

    def generic(self, args, kws):
        cont, item = args
        if isinstance(cont, types.Container):
            return signature(types.boolean, cont, cont.dtype)

@infer_global(len)
class ContainerLen(AbstractTemplate):

    def generic(self, args, kws):
        assert not kws
        (val,) = args
        if isinstance(val, (types.Container)):
            return signature(types.intp, val)


@infer_global(operator.truth)
class SequenceBool(AbstractTemplate):
    key = operator.truth

    def generic(self, args, kws):
        assert not kws
        (val,) = args
        if isinstance(val, (types.Sequence)):
            return signature(types.boolean, val)


@infer_global(operator.getitem)
class GetItemSequence(AbstractTemplate):
    key = operator.getitem

    def generic(self, args, kws):
        seq, idx = args
        if isinstance(seq, types.Sequence):
            idx = normalize_1d_index(idx)
            if isinstance(idx, types.SliceType):
                # Slicing a tuple only supported with static_getitem
                if not isinstance(seq, types.BaseTuple):
                    return signature(seq, seq, idx)
            elif isinstance(idx, types.Integer):
                return signature(seq.dtype, seq, idx)

@infer_global(operator.setitem)
class SetItemSequence(AbstractTemplate):
    def generic(self, args, kws):
        seq, idx, value = args
        if isinstance(seq, types.MutableSequence):
            idx = normalize_1d_index(idx)
            if isinstance(idx, types.SliceType):
                return signature(types.none, seq, idx, seq)
            elif isinstance(idx, types.Integer):
                if not self.context.can_convert(value, seq.dtype):
                    msg = "invalid setitem with value of {} to element of {}"
                    raise errors.TypingError(msg.format(types.unliteral(value), seq.dtype))
                return signature(types.none, seq, idx, seq.dtype)


@infer_global(operator.delitem)
class DelItemSequence(AbstractTemplate):
    def generic(self, args, kws):
        seq, idx = args
        if isinstance(seq, types.MutableSequence):
            idx = normalize_1d_index(idx)
            return signature(types.none, seq, idx)


# --------------------------------------------------------------------------
# named tuples

@infer_getattr
class NamedTupleAttribute(AttributeTemplate):
    key = types.BaseNamedTuple

    def resolve___class__(self, tup):
        return types.NamedTupleClass(tup.instance_class)

    def generic_resolve(self, tup, attr):
        # Resolution of other attributes
        try:
            index = tup.fields.index(attr)
        except ValueError:
            return
        return tup[index]


@infer_getattr
class NamedTupleClassAttribute(AttributeTemplate):
    key = types.NamedTupleClass

    def resolve___call__(self, classty):
        """
        Resolve the named tuple constructor, aka the class's __call__ method.
        """
        instance_class = classty.instance_class
        pysig = utils.pysignature(instance_class)

        def typer(*args, **kws):
            # Fold keyword args
            try:
                bound = pysig.bind(*args, **kws)
            except TypeError as e:
                msg = "In '%s': %s" % (instance_class, e)
                e.args = (msg,)
                raise
            assert not bound.kwargs
            return types.BaseTuple.from_types(bound.args, instance_class)

        # Override the typer's pysig to match the namedtuple constructor's
        typer.pysig = pysig
        return types.Function(make_callable_template(self.key, typer))
