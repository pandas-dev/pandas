from numba.core import errors, ir, consts
from numba.core.rewrites import register_rewrite, Rewrite


@register_rewrite('before-inference')
class RewriteConstRaises(Rewrite):
    """
    Rewrite IR statements of the kind `raise(value)`
    where `value` is the result of instantiating an exception with
    constant arguments
    into `static_raise(exception_type, constant args)`.

    This allows lowering in nopython mode, where one can't instantiate
    exception instances from runtime data.
    """

    def _is_exception_type(self, const):
        return isinstance(const, type) and issubclass(const, Exception)

    def _break_constant(self, const, loc):
        """
        Break down constant exception.
        """
        if isinstance(const, tuple): # it's a tuple(exception class, args)
            if not self._is_exception_type(const[0]):
                msg = "Encountered unsupported exception constant %r"
                raise errors.UnsupportedError(msg % (const[0],), loc)
            return const[0], tuple(const[1])
        elif self._is_exception_type(const):
            return const, None
        else:
            if isinstance(const, str):
                msg = ("Directly raising a string constant as an exception is "
                       "not supported.")
            else:
                msg = "Encountered unsupported constant type used for exception"
            raise errors.UnsupportedError(msg, loc)

    def _try_infer_constant(self, func_ir, inst):
        try:
            return func_ir.infer_constant(inst.exception)
        except consts.ConstantInferenceError:
            # not a static exception
            return None

    def match(self, func_ir, block, typemap, calltypes):
        self.raises = raises = {}
        self.tryraises = tryraises = {}
        self.block = block
        # Detect all raise statements and find which ones can be
        # rewritten
        for inst in block.find_insts((ir.Raise, ir.TryRaise)):
            if inst.exception is None:
                # re-reraise
                exc_type, exc_args = None, None
            else:
                # raise <something> => find the definition site for <something>
                const = self._try_infer_constant(func_ir, inst)

                # failure to infer constant indicates this isn't a static
                # exception
                if const is None:
                    continue

                loc = inst.exception.loc
                exc_type, exc_args = self._break_constant(const, loc)

            if isinstance(inst, ir.Raise):
                raises[inst] = exc_type, exc_args
            elif isinstance(inst, ir.TryRaise):
                tryraises[inst] = exc_type, exc_args
            else:
                raise ValueError('unexpected: {}'.format(type(inst)))
        return (len(raises) + len(tryraises)) > 0

    def apply(self):
        """
        Rewrite all matching setitems as static_setitems.
        """
        new_block = self.block.copy()
        new_block.clear()
        for inst in self.block.body:
            if inst in self.raises:
                exc_type, exc_args = self.raises[inst]
                new_inst = ir.StaticRaise(exc_type, exc_args, inst.loc)
                new_block.append(new_inst)
            elif inst in self.tryraises:
                exc_type, exc_args = self.tryraises[inst]
                new_inst = ir.StaticTryRaise(exc_type, exc_args, inst.loc)
                new_block.append(new_inst)
            else:
                new_block.append(inst)
        return new_block
