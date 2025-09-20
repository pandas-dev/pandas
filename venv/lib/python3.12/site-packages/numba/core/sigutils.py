from numba.core import types, typing


def is_signature(sig):
    """
    Return whether *sig* is a potentially valid signature
    specification (for user-facing APIs).
    """
    return isinstance(sig, (str, tuple, typing.Signature))


def _parse_signature_string(signature_str):
    """
    Parameters
    ----------
    signature_str : str
    """
    # Just eval signature_str using the types submodules as globals
    return eval(signature_str, {}, types.__dict__)


def normalize_signature(sig):
    """
    From *sig* (a signature specification), return a ``(args, return_type)``
    tuple, where ``args`` itself is a tuple of types, and ``return_type``
    can be None if not specified.
    """
    if isinstance(sig, str):
        parsed = _parse_signature_string(sig)
    else:
        parsed = sig
    if isinstance(parsed, tuple):
        args, return_type = parsed, None
    elif isinstance(parsed, typing.Signature):
        args, return_type = parsed.args, parsed.return_type
    else:
        raise TypeError("invalid signature: %r (type: %r) evaluates to %r "
                        "instead of tuple or Signature" % (
                            sig, sig.__class__.__name__,
                            parsed.__class__.__name__
                        ))

    def check_type(ty):
        if not isinstance(ty, types.Type):
            raise TypeError("invalid type in signature: expected a type "
                            "instance, got %r" % (ty,))

    if return_type is not None:
        check_type(return_type)
    for ty in args:
        check_type(ty)

    return args, return_type
