
__all__ = ['FunctionType', 'UndefinedFunctionType', 'FunctionPrototype',
           'WrapperAddressProtocol', 'CompileResultWAP']

from abc import ABC, abstractmethod
from .abstract import Type
from .. import types, errors


class FunctionType(Type):
    """
    First-class function type.
    """

    cconv = None

    def __init__(self, signature):
        sig = types.unliteral(signature)
        self.nargs = len(sig.args)
        self.signature = sig
        self.ftype = FunctionPrototype(sig.return_type, sig.args)
        self._key = self.ftype.key

    @property
    def key(self):
        return self._key

    @property
    def name(self):
        return f'{type(self).__name__}[{self.key}]'

    def is_precise(self):
        return self.signature.is_precise()

    def get_precise(self):
        return self

    def dump(self, tab=''):
        print(f'{tab}DUMP {type(self).__name__}[code={self._code}]')
        self.signature.dump(tab=tab + '  ')
        print(f'{tab}END DUMP {type(self).__name__}')

    def get_call_type(self, context, args, kws):
        from numba.core import typing

        if kws:
            # First-class functions carry only the type signature
            # information and function address value. So, it is not
            # possible to determine the positional arguments
            # corresponding to the keyword arguments in the call
            # expression. For instance, the definition of the
            # first-class function may not use the same argument names
            # that the caller assumes. [numba/issues/5540].
            raise errors.UnsupportedError(
                'first-class function call cannot use keyword arguments')

        if len(args) != self.nargs:
            raise ValueError(
                f'mismatch of arguments number: {len(args)} vs {self.nargs}')

        sig = self.signature

        # check that arguments types match with the signature types exactly
        for atype, sig_atype in zip(args, sig.args):
            atype = types.unliteral(atype)
            if sig_atype.is_precise():
                conv_score = context.context.can_convert(
                    fromty=atype, toty=sig_atype
                )
                if conv_score is None \
                   or conv_score > typing.context.Conversion.safe:
                    raise ValueError(
                        f'mismatch of argument types: {atype} vs {sig_atype}')

        if not sig.is_precise():
            for dispatcher in self.dispatchers:
                template, pysig, args, kws \
                    = dispatcher.get_call_template(args, kws)
                new_sig = template(context.context).apply(args, kws)
                return types.unliteral(new_sig)

        return sig

    def check_signature(self, other_sig):
        """Return True if signatures match (up to being precise).
        """
        sig = self.signature
        return (self.nargs == len(other_sig.args)
                and (sig == other_sig or not sig.is_precise()))

    def unify(self, context, other):
        if isinstance(other, types.UndefinedFunctionType) \
           and self.nargs == other.nargs:
            return self


class UndefinedFunctionType(FunctionType):

    _counter = 0

    def __init__(self, nargs, dispatchers):
        from numba.core.typing.templates import Signature
        signature = Signature(types.undefined,
                              (types.undefined,) * nargs, recvr=None)

        super(UndefinedFunctionType, self).__init__(signature)

        self.dispatchers = dispatchers

        # make the undefined function type instance unique
        type(self)._counter += 1
        self._key += str(type(self)._counter)

    def get_precise(self):
        """
        Return precise function type if possible.
        """
        for dispatcher in self.dispatchers:
            for cres in dispatcher.overloads.values():
                sig = types.unliteral(cres.signature)
                return FunctionType(sig)
        return self


class FunctionPrototype(Type):
    """
    Represents the prototype of a first-class function type.
    Used internally.
    """
    cconv = None

    def __init__(self, rtype, atypes):
        self.rtype = rtype
        self.atypes = tuple(atypes)

        assert isinstance(rtype, Type), (rtype)
        lst = []
        for atype in self.atypes:
            assert isinstance(atype, Type), (atype)
            lst.append(atype.name)
        name = '%s(%s)' % (rtype, ', '.join(lst))

        super(FunctionPrototype, self).__init__(name)

    @property
    def key(self):
        return self.name


class WrapperAddressProtocol(ABC):
    """Base class for Wrapper Address Protocol.

    Objects that inherit from the WrapperAddressProtocol can be passed
    as arguments to Numba jit compiled functions where it can be used
    as first-class functions. As a minimum, the derived types must
    implement two methods ``__wrapper_address__`` and ``signature``.
    """

    @abstractmethod
    def __wrapper_address__(self):
        """Return the address of a first-class function.

        Returns
        -------
        addr : int
        """

    @abstractmethod
    def signature(self):
        """Return the signature of a first-class function.

        Returns
        -------
        sig : Signature
          The returned Signature instance represents the type of a
          first-class function that the given WrapperAddressProtocol
          instance represents.
        """


class CompileResultWAP(WrapperAddressProtocol):
    """Wrapper of dispatcher instance compilation result to turn it a
    first-class function.
    """

    def __init__(self, cres):
        """
        Parameters
        ----------
        cres : CompileResult
          Specify compilation result of a Numba jit-decorated function
          (that is a value of dispatcher instance ``overloads``
          attribute)
        """
        self.cres = cres
        name = getattr(cres.fndesc, 'llvm_cfunc_wrapper_name')
        self.address = cres.library.get_pointer_to_function(name)

    def dump(self, tab=''):
        print(f'{tab}DUMP {type(self).__name__} [addr={self.address}]')
        self.cres.signature.dump(tab=tab + '  ')
        print(f'{tab}END DUMP {type(self).__name__}')

    def __wrapper_address__(self):
        return self.address

    def signature(self):
        return self.cres.signature

    def __call__(self, *args, **kwargs):  # used in object-mode
        return self.cres.entry_point(*args, **kwargs)
