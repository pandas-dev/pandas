"""
Define typing templates
"""

from abc import ABC, abstractmethod
import functools
import sys
import inspect
import os.path
from collections import namedtuple
from collections.abc import Sequence
from types import MethodType, FunctionType, MappingProxyType

import numba
from numba.core import types, utils, targetconfig
from numba.core.errors import (
    TypingError,
    InternalError,
)
from numba.core.cpu_options import InlineOptions

# info store for inliner callback functions e.g. cost model
_inline_info = namedtuple('inline_info',
                          'func_ir typemap calltypes signature')


class Signature(object):
    """
    The signature of a function call or operation, i.e. its argument types
    and return type.
    """

    # XXX Perhaps the signature should be a BoundArguments, instead
    # of separate args and pysig...
    __slots__ = '_return_type', '_args', '_recvr', '_pysig'

    def __init__(self, return_type, args, recvr, pysig=None):
        if isinstance(args, list):
            args = tuple(args)
        self._return_type = return_type
        self._args = args
        self._recvr = recvr
        self._pysig = pysig

    @property
    def return_type(self):
        return self._return_type

    @property
    def args(self):
        return self._args

    @property
    def recvr(self):
        return self._recvr

    @property
    def pysig(self):
        return self._pysig

    def replace(self, **kwargs):
        """Copy and replace the given attributes provided as keyword arguments.
        Returns an updated copy.
        """
        curstate = dict(return_type=self.return_type,
                        args=self.args,
                        recvr=self.recvr,
                        pysig=self.pysig)
        curstate.update(kwargs)
        return Signature(**curstate)

    def __getstate__(self):
        """
        Needed because of __slots__.
        """
        return self._return_type, self._args, self._recvr, self._pysig

    def __setstate__(self, state):
        """
        Needed because of __slots__.
        """
        self._return_type, self._args, self._recvr, self._pysig = state

    def __hash__(self):
        return hash((self.args, self.return_type))

    def __eq__(self, other):
        if isinstance(other, Signature):
            return (self.args == other.args and
                    self.return_type == other.return_type and
                    self.recvr == other.recvr and
                    self.pysig == other.pysig)

    def __ne__(self, other):
        return not (self == other)

    def __repr__(self):
        return "%s -> %s" % (self.args, self.return_type)

    @property
    def is_method(self):
        """
        Whether this signature represents a bound method or a regular
        function.
        """
        return self.recvr is not None

    def as_method(self):
        """
        Convert this signature to a bound method signature.
        """
        if self.recvr is not None:
            return self
        sig = signature(self.return_type, *self.args[1:],
                        recvr=self.args[0])

        # Adjust the python signature
        params = list(self.pysig.parameters.values())[1:]
        sig = sig.replace(
            pysig=utils.pySignature(
                parameters=params,
                return_annotation=self.pysig.return_annotation,
            ),
        )
        return sig

    def as_function(self):
        """
        Convert this signature to a regular function signature.
        """
        if self.recvr is None:
            return self
        sig = signature(self.return_type, *((self.recvr,) + self.args))
        return sig

    def as_type(self):
        """
        Convert this signature to a first-class function type.
        """
        return types.FunctionType(self)

    def __unliteral__(self):
        return signature(types.unliteral(self.return_type),
                         *map(types.unliteral, self.args))

    def dump(self, tab=''):
        c = self.as_type()._code
        print(f'{tab}DUMP {type(self).__name__} [type code: {c}]')
        print(f'{tab}  Argument types:')
        for a in self.args:
            a.dump(tab=tab + '  | ')
        print(f'{tab}  Return type:')
        self.return_type.dump(tab=tab + '  | ')
        print(f'{tab}END DUMP')

    def is_precise(self):
        for atype in self.args:
            if not atype.is_precise():
                return False
        return self.return_type.is_precise()


def make_concrete_template(name, key, signatures):
    baseclasses = (ConcreteTemplate,)
    gvars = dict(key=key, cases=list(signatures))
    return type(name, baseclasses, gvars)


def make_callable_template(key, typer, recvr=None):
    """
    Create a callable template with the given key and typer function.
    """
    def generic(self):
        return typer

    name = "%s_CallableTemplate" % (key,)
    bases = (CallableTemplate,)
    class_dict = dict(key=key, generic=generic, recvr=recvr)
    return type(name, bases, class_dict)


def signature(return_type, *args, **kws):
    recvr = kws.pop('recvr', None)
    assert not kws
    return Signature(return_type, args, recvr=recvr)


def fold_arguments(pysig, args, kws, normal_handler, default_handler,
                   stararg_handler):
    """
    Given the signature *pysig*, explicit *args* and *kws*, resolve
    omitted arguments and keyword arguments. A tuple of positional
    arguments is returned.
    Various handlers allow to process arguments:
    - normal_handler(index, param, value) is called for normal arguments
    - default_handler(index, param, default) is called for omitted arguments
    - stararg_handler(index, param, values) is called for a "*args" argument
    """
    if isinstance(kws, Sequence):
        # Normalize dict kws
        kws = dict(kws)

    # deal with kwonly args
    params = pysig.parameters
    kwonly = []
    for name, p in params.items():
        if p.kind == p.KEYWORD_ONLY:
            kwonly.append(name)

    if kwonly:
        bind_args = args[:-len(kwonly)]
    else:
        bind_args = args
    bind_kws = kws.copy()
    if kwonly:
        for idx, n in enumerate(kwonly):
            bind_kws[n] = args[len(kwonly) + idx]

    # now bind
    ba = pysig.bind(*bind_args, **bind_kws)
    for i, param in enumerate(pysig.parameters.values()):
        name = param.name
        default = param.default
        if param.kind == param.VAR_POSITIONAL:
            # stararg may be omitted, in which case its "default" value
            # is simply the empty tuple
            if name in ba.arguments:
                argval = ba.arguments[name]
                # NOTE: avoid wrapping the tuple type for stararg in another
                #       tuple.
                if (len(argval) == 1 and
                        isinstance(argval[0], (types.StarArgTuple,
                                               types.StarArgUniTuple))):
                    argval = tuple(argval[0])
            else:
                argval = ()
            out = stararg_handler(i, param, argval)

            ba.arguments[name] = out
        elif name in ba.arguments:
            # Non-stararg, present
            ba.arguments[name] = normal_handler(i, param, ba.arguments[name])
        else:
            # Non-stararg, omitted
            assert default is not param.empty
            ba.arguments[name] = default_handler(i, param, default)
    # Collect args in the right order
    args = tuple(ba.arguments[param.name]
                 for param in pysig.parameters.values())
    return args


class FunctionTemplate(ABC):
    # Set to true to disable unsafe cast.
    # subclass overide-able
    unsafe_casting = True
    # Set to true to require exact match without casting.
    # subclass overide-able
    exact_match_required = False
    # Set to true to prefer literal arguments.
    # Useful for definitions that specialize on literal but also support
    # non-literals.
    # subclass overide-able
    prefer_literal = False
    # metadata
    metadata = {}

    def __init__(self, context):
        self.context = context

    def _select(self, cases, args, kws):
        options = {
            'unsafe_casting': self.unsafe_casting,
            'exact_match_required': self.exact_match_required,
        }
        selected = self.context.resolve_overload(self.key, cases, args, kws,
                                                 **options)
        return selected

    def get_impl_key(self, sig):
        """
        Return the key for looking up the implementation for the given
        signature on the target context.
        """
        # Lookup the key on the class, to avoid binding it with `self`.
        key = type(self).key
        # On Python 2, we must also take care about unbound methods
        if isinstance(key, MethodType):
            assert key.im_self is None
            key = key.im_func
        return key

    @classmethod
    def get_source_code_info(cls, impl):
        """
        Gets the source information about function impl.
        Returns:

        code - str: source code as a string
        firstlineno - int: the first line number of the function impl
        path - str: the path to file containing impl

        if any of the above are not available something generic is returned
        """
        try:
            code, firstlineno = inspect.getsourcelines(impl)
        except OSError: # missing source, probably a string
            code = "None available (built from string?)"
            firstlineno = 0
        path = inspect.getsourcefile(impl)
        if path is None:
            path = "<unknown> (built from string?)"
        return code, firstlineno, path

    @abstractmethod
    def get_template_info(self):
        """
        Returns a dictionary with information specific to the template that will
        govern how error messages are displayed to users. The dictionary must
        be of the form:
        info = {
            'kind': "unknown", # str: The kind of template, e.g. "Overload"
            'name': "unknown", # str: The name of the source function
            'sig': "unknown",  # str: The signature(s) of the source function
            'filename': "unknown", # str: The filename of the source function
            'lines': ("start", "end"), # tuple(int, int): The start and
                                         end line of the source function.
            'docstring': "unknown" # str: The docstring of the source function
        }
        """
        pass

    def __str__(self):
        info = self.get_template_info()
        srcinfo = f"{info['filename']}:{info['lines'][0]}"
        return f"<{self.__class__.__name__} {srcinfo}>"

    __repr__ = __str__


class AbstractTemplate(FunctionTemplate):
    """
    Defines method ``generic(self, args, kws)`` which compute a possible
    signature base on input types.  The signature does not have to match the
    input types. It is compared against the input types afterwards.
    """

    def apply(self, args, kws):
        generic = getattr(self, "generic")
        sig = generic(args, kws)
        # Enforce that *generic()* must return None or Signature
        if sig is not None:
            if not isinstance(sig, Signature):
                raise AssertionError(
                    "generic() must return a Signature or None. "
                    "{} returned {}".format(generic, type(sig)),
                )

        # Unpack optional type if no matching signature
        if not sig and any(isinstance(x, types.Optional) for x in args):
            def unpack_opt(x):
                if isinstance(x, types.Optional):
                    return x.type
                else:
                    return x

            args = list(map(unpack_opt, args))
            assert not kws  # Not supported yet
            sig = generic(args, kws)

        return sig

    def get_template_info(self):
        impl = getattr(self, "generic")
        basepath = os.path.dirname(os.path.dirname(numba.__file__))

        code, firstlineno, path = self.get_source_code_info(impl)
        sig = str(utils.pysignature(impl))
        info = {
            'kind': "overload",
            'name': getattr(impl, '__qualname__', impl.__name__),
            'sig': sig,
            'filename': utils.safe_relpath(path, start=basepath),
            'lines': (firstlineno, firstlineno + len(code) - 1),
            'docstring': impl.__doc__
        }
        return info


class CallableTemplate(FunctionTemplate):
    """
    Base class for a template defining a ``generic(self)`` method
    returning a callable to be called with the actual ``*args`` and
    ``**kwargs`` representing the call signature.  The callable has
    to return a return type, a full signature, or None.  The signature
    does not have to match the input types. It is compared against the
    input types afterwards.
    """
    recvr = None

    def apply(self, args, kws):
        generic = getattr(self, "generic")
        typer = generic()
        match_sig = inspect.signature(typer)
        try:
            match_sig.bind(*args, **kws)
        except TypeError as e:
            # bind failed, raise, if there's a
            # ValueError then there's likely unrecoverable
            # problems
            raise TypingError(str(e)) from e

        sig = typer(*args, **kws)

        # Unpack optional type if no matching signature
        if sig is None:
            if any(isinstance(x, types.Optional) for x in args):
                def unpack_opt(x):
                    if isinstance(x, types.Optional):
                        return x.type
                    else:
                        return x

                args = list(map(unpack_opt, args))
                sig = typer(*args, **kws)
            if sig is None:
                return

        # Get the pysig
        try:
            pysig = typer.pysig
        except AttributeError:
            pysig = utils.pysignature(typer)

        # Fold any keyword arguments
        bound = pysig.bind(*args, **kws)
        if bound.kwargs:
            raise TypingError("unsupported call signature")
        if not isinstance(sig, Signature):
            # If not a signature, `sig` is assumed to be the return type
            if not isinstance(sig, types.Type):
                raise TypeError("invalid return type for callable template: "
                                "got %r" % (sig,))
            sig = signature(sig, *bound.args)
        if self.recvr is not None:
            sig = sig.replace(recvr=self.recvr)
        # Hack any omitted parameters out of the typer's pysig,
        # as lowering expects an exact match between formal signature
        # and actual args.
        if len(bound.args) < len(pysig.parameters):
            parameters = list(pysig.parameters.values())[:len(bound.args)]
            pysig = pysig.replace(parameters=parameters)
        sig = sig.replace(pysig=pysig)
        cases = [sig]
        return self._select(cases, bound.args, bound.kwargs)

    def get_template_info(self):
        impl = getattr(self, "generic")
        basepath = os.path.dirname(os.path.dirname(numba.__file__))
        code, firstlineno, path = self.get_source_code_info(impl)
        sig = str(utils.pysignature(impl))
        info = {
            'kind': "overload",
            'name': getattr(self.key, '__name__',
                            getattr(impl, '__qualname__', impl.__name__),),
            'sig': sig,
            'filename': utils.safe_relpath(path, start=basepath),
            'lines': (firstlineno, firstlineno + len(code) - 1),
            'docstring': impl.__doc__
        }
        return info


class ConcreteTemplate(FunctionTemplate):
    """
    Defines attributes "cases" as a list of signature to match against the
    given input types.
    """

    def apply(self, args, kws):
        cases = getattr(self, 'cases')
        return self._select(cases, args, kws)

    def get_template_info(self):
        import operator
        name = getattr(self.key, '__name__', "unknown")
        op_func = getattr(operator, name, None)

        kind = "Type restricted function"
        if op_func is not None:
            if self.key is op_func:
                kind = "operator overload"
        info = {
            'kind': kind,
            'name': name,
            'sig': "unknown",
            'filename': "unknown",
            'lines': ("unknown", "unknown"),
            'docstring': "unknown"
        }
        return info


class _EmptyImplementationEntry(InternalError):
    def __init__(self, reason):
        super(_EmptyImplementationEntry, self).__init__(
            "_EmptyImplementationEntry({!r})".format(reason),
        )


class _OverloadFunctionTemplate(AbstractTemplate):
    """
    A base class of templates for overload functions.
    """

    def _validate_sigs(self, typing_func, impl_func):
        # check that the impl func and the typing func have the same signature!
        typing_sig = utils.pysignature(typing_func)
        impl_sig = utils.pysignature(impl_func)
        # the typing signature is considered golden and must be adhered to by
        # the implementation...
        # Things that are valid:
        # 1. args match exactly
        # 2. kwargs match exactly in name and default value
        # 3. Use of *args in the same location by the same name in both typing
        #    and implementation signature
        # 4. Use of *args in the implementation signature to consume any number
        #    of arguments in the typing signature.
        # Things that are invalid:
        # 5. Use of *args in the typing signature that is not replicated
        #    in the implementing signature
        # 6. Use of **kwargs

        def get_args_kwargs(sig):
            kws = []
            args = []
            pos_arg = None
            for x in sig.parameters.values():
                if x.default == utils.pyParameter.empty:
                    args.append(x)
                    if x.kind == utils.pyParameter.VAR_POSITIONAL:
                        pos_arg = x
                    elif x.kind == utils.pyParameter.VAR_KEYWORD:
                        msg = ("The use of VAR_KEYWORD (e.g. **kwargs) is "
                               "unsupported. (offending argument name is '%s')")
                        raise InternalError(msg % x)
                else:
                    kws.append(x)
            return args, kws, pos_arg

        ty_args, ty_kws, ty_pos = get_args_kwargs(typing_sig)
        im_args, im_kws, im_pos = get_args_kwargs(impl_sig)

        sig_fmt = ("Typing signature:         %s\n"
                   "Implementation signature: %s")
        sig_str = sig_fmt % (typing_sig, impl_sig)

        err_prefix = "Typing and implementation arguments differ in "

        a = ty_args
        b = im_args
        if ty_pos:
            if not im_pos:
                # case 5. described above
                msg = ("VAR_POSITIONAL (e.g. *args) argument kind (offending "
                       "argument name is '%s') found in the typing function "
                       "signature, but is not in the implementing function "
                       "signature.\n%s") % (ty_pos, sig_str)
                raise InternalError(msg)
        else:
            if im_pos:
                # no *args in typing but there's a *args in the implementation
                # this is case 4. described above
                b = im_args[:im_args.index(im_pos)]
                try:
                    a = ty_args[:ty_args.index(b[-1]) + 1]
                except ValueError:
                    # there's no b[-1] arg name in the ty_args, something is
                    # very wrong, we can't work out a diff (*args consumes
                    # unknown quantity of args) so just report first error
                    specialized = "argument names.\n%s\nFirst difference: '%s'"
                    msg = err_prefix + specialized % (sig_str, b[-1])
                    raise InternalError(msg)

        def gen_diff(typing, implementing):
            diff = set(typing) ^ set(implementing)
            return "Difference: %s" % diff

        if a != b:
            specialized = "argument names.\n%s\n%s" % (sig_str, gen_diff(a, b))
            raise InternalError(err_prefix + specialized)

        # ensure kwargs are the same
        ty = [x.name for x in ty_kws]
        im = [x.name for x in im_kws]
        if ty != im:
            specialized = "keyword argument names.\n%s\n%s"
            msg = err_prefix + specialized % (sig_str, gen_diff(ty_kws, im_kws))
            raise InternalError(msg)
        same = [x.default for x in ty_kws] == [x.default for x in im_kws]
        if not same:
            specialized = "keyword argument default values.\n%s\n%s"
            msg = err_prefix + specialized % (sig_str, gen_diff(ty_kws, im_kws))
            raise InternalError(msg)

    def generic(self, args, kws):
        """
        Type the overloaded function by compiling the appropriate
        implementation for the given args.
        """
        from numba.core.typed_passes import PreLowerStripPhis

        disp, new_args = self._get_impl(args, kws)
        if disp is None:
            return
        # Compile and type it for the given types
        disp_type = types.Dispatcher(disp)
        # Store the compiled overload for use in the lowering phase if there's
        # no inlining required (else functions are being compiled which will
        # never be used as they are inlined)
        if not self._inline.is_never_inline:
            # need to run the compiler front end up to type inference to compute
            # a signature
            from numba.core import typed_passes, compiler
            from numba.core.inline_closurecall import InlineWorker
            fcomp = disp._compiler
            flags = compiler.Flags()

            # Updating these causes problems?!
            #fcomp.targetdescr.options.parse_as_flags(flags,
            #                                         fcomp.targetoptions)
            #flags = fcomp._customize_flags(flags)

            # spoof a compiler pipline like the one that will be in use
            tyctx = fcomp.targetdescr.typing_context
            tgctx = fcomp.targetdescr.target_context
            compiler_inst = fcomp.pipeline_class(tyctx, tgctx, None, None, None,
                                                 flags, None, )
            inline_worker = InlineWorker(tyctx, tgctx, fcomp.locals,
                                         compiler_inst, flags, None,)

            # If the inlinee contains something to trigger literal arg dispatch
            # then the pipeline call will unconditionally fail due to a raised
            # ForceLiteralArg exception. Therefore `resolve` is run first, as
            # type resolution must occur at some point, this will hit any
            # `literally` calls and because it's going via the dispatcher will
            # handle them correctly i.e. ForceLiteralArg propagates. This having
            # the desired effect of ensuring the pipeline call is only made in
            # situations that will succeed. For context see #5887.
            resolve = disp_type.dispatcher.get_call_template
            template, pysig, folded_args, kws = resolve(new_args, kws)
            ir = inline_worker.run_untyped_passes(
                disp_type.dispatcher.py_func, enable_ssa=True
            )

            (
                typemap,
                return_type,
                calltypes,
                _
            ) = typed_passes.type_inference_stage(
                self.context, tgctx, ir, folded_args, None)
            ir = PreLowerStripPhis()._strip_phi_nodes(ir)
            ir._definitions = numba.core.ir_utils.build_definitions(ir.blocks)

            sig = Signature(return_type, folded_args, None)
            # this stores a load of info for the cost model function if supplied
            # it by default is None
            self._inline_overloads[sig.args] = {'folded_args': folded_args}
            # this stores the compiled overloads, if there's no compiled
            # overload available i.e. function is always inlined, the key still
            # needs to exist for type resolution

            # NOTE: If lowering is failing on a `_EmptyImplementationEntry`,
            #       the inliner has failed to inline this entry correctly.
            impl_init = _EmptyImplementationEntry('always inlined')
            self._compiled_overloads[sig.args] = impl_init
            if not self._inline.is_always_inline:
                # this branch is here because a user has supplied a function to
                # determine whether to inline or not. As a result both compiled
                # function and inliner info needed, delaying the computation of
                # this leads to an internal state mess at present. TODO: Fix!
                sig = disp_type.get_call_type(self.context, new_args, kws)
                self._compiled_overloads[sig.args] = disp_type.get_overload(sig)
                # store the inliner information, it's used later in the cost
                # model function call
            iinfo = _inline_info(ir, typemap, calltypes, sig)
            self._inline_overloads[sig.args] = {'folded_args': folded_args,
                                                'iinfo': iinfo}
        else:
            sig = disp_type.get_call_type(self.context, new_args, kws)
            if sig is None: # can't resolve for this target
                return None
            self._compiled_overloads[sig.args] = disp_type.get_overload(sig)
        return sig

    def _get_impl(self, args, kws):
        """Get implementation given the argument types.

        Returning a Dispatcher object.  The Dispatcher object is cached
        internally in `self._impl_cache`.
        """
        flags = targetconfig.ConfigStack.top_or_none()
        cache_key = self.context, tuple(args), tuple(kws.items()), flags
        try:
            impl, args = self._impl_cache[cache_key]
            return impl, args
        except KeyError:
            # pass and try outside the scope so as to not have KeyError with a
            # nested addition error in the case the _build_impl fails
            pass
        impl, args = self._build_impl(cache_key, args, kws)
        return impl, args

    def _get_jit_decorator(self):
        """Gets a jit decorator suitable for the current target"""

        from numba.core.target_extension import (target_registry,
                                                 get_local_target,
                                                 jit_registry)

        jitter_str = self.metadata.get('target', 'generic')
        jitter = jit_registry.get(jitter_str, None)

        if jitter is None:
            # No JIT known for target string, see if something is
            # registered for the string and report if not.
            target_class = target_registry.get(jitter_str, None)
            if target_class is None:
                msg = ("Unknown target '{}', has it been ",
                       "registered?")
                raise ValueError(msg.format(jitter_str))

            target_hw = get_local_target(self.context)

            # check that the requested target is in the hierarchy for the
            # current frame's target.
            if not issubclass(target_hw, target_class):
                msg = "No overloads exist for the requested target: {}."

            jitter = jit_registry[target_hw]

        if jitter is None:
            raise ValueError("Cannot find a suitable jit decorator")

        return jitter

    def _build_impl(self, cache_key, args, kws):
        """Build and cache the implementation.

        Given the positional (`args`) and keyword arguments (`kws`), obtains
        the `overload` implementation and wrap it in a Dispatcher object.
        The expected argument types are returned for use by type-inference.
        The expected argument types are only different from the given argument
        types if there is an imprecise type in the given argument types.

        Parameters
        ----------
        cache_key : hashable
            The key used for caching the implementation.
        args : Tuple[Type]
            Types of positional argument.
        kws : Dict[Type]
            Types of keyword argument.

        Returns
        -------
        disp, args :
            On success, returns `(Dispatcher, Tuple[Type])`.
            On failure, returns `(None, None)`.

        """
        jitter = self._get_jit_decorator()

        # Get the overload implementation for the given types
        ov_sig = inspect.signature(self._overload_func)
        try:
            ov_sig.bind(*args, **kws)
        except TypeError as e:
            # bind failed, raise, if there's a
            # ValueError then there's likely unrecoverable
            # problems
            raise TypingError(str(e)) from e
        else:
            ovf_result = self._overload_func(*args, **kws)

        if ovf_result is None:
            # No implementation => fail typing
            self._impl_cache[cache_key] = None, None
            return None, None
        elif isinstance(ovf_result, tuple):
            # The implementation returned a signature that the type-inferencer
            # should be using.
            sig, pyfunc = ovf_result
            args = sig.args
            kws = {}
            cache_key = None            # don't cache
        else:
            # Regular case
            pyfunc = ovf_result

        # Check type of pyfunc
        if not isinstance(pyfunc, FunctionType):
            msg = ("Implementation function returned by `@overload` "
                   "has an unexpected type.  Got {}")
            raise AssertionError(msg.format(pyfunc))

        # check that the typing and impl sigs match up
        if self._strict:
            self._validate_sigs(self._overload_func, pyfunc)
        # Make dispatcher
        jitdecor = jitter(**self._jit_options)
        disp = jitdecor(pyfunc)
        # Make sure that the implementation can be fully compiled
        disp_type = types.Dispatcher(disp)
        disp_type.get_call_type(self.context, args, kws)
        if cache_key is not None:
            self._impl_cache[cache_key] = disp, args
        return disp, args

    def get_impl_key(self, sig):
        """
        Return the key for looking up the implementation for the given
        signature on the target context.
        """
        return self._compiled_overloads[sig.args]

    @classmethod
    def get_source_info(cls):
        """Return a dictionary with information about the source code of the
        implementation.

        Returns
        -------
        info : dict
            - "kind" : str
                The implementation kind.
            - "name" : str
                The name of the function that provided the definition.
            - "sig" : str
                The formatted signature of the function.
            - "filename" : str
                The name of the source file.
            - "lines": tuple (int, int)
                First and list line number.
            - "docstring": str
                The docstring of the definition.
        """
        basepath = os.path.dirname(os.path.dirname(numba.__file__))
        impl = cls._overload_func
        code, firstlineno, path = cls.get_source_code_info(impl)
        sig = str(utils.pysignature(impl))
        info = {
            'kind': "overload",
            'name': getattr(impl, '__qualname__', impl.__name__),
            'sig': sig,
            'filename': utils.safe_relpath(path, start=basepath),
            'lines': (firstlineno, firstlineno + len(code) - 1),
            'docstring': impl.__doc__
        }
        return info

    def get_template_info(self):
        basepath = os.path.dirname(os.path.dirname(numba.__file__))
        impl = self._overload_func
        code, firstlineno, path = self.get_source_code_info(impl)
        sig = str(utils.pysignature(impl))
        info = {
            'kind': "overload",
            'name': getattr(impl, '__qualname__', impl.__name__),
            'sig': sig,
            'filename': utils.safe_relpath(path, start=basepath),
            'lines': (firstlineno, firstlineno + len(code) - 1),
            'docstring': impl.__doc__
        }
        return info


def make_overload_template(func, overload_func, jit_options, strict,
                           inline, prefer_literal=False, **kwargs):
    """
    Make a template class for function *func* overloaded by *overload_func*.
    Compiler options are passed as a dictionary to *jit_options*.
    """
    func_name = getattr(func, '__name__', str(func))
    name = "OverloadTemplate_%s" % (func_name,)
    base = _OverloadFunctionTemplate
    dct = dict(key=func, _overload_func=staticmethod(overload_func),
               _impl_cache={}, _compiled_overloads={}, _jit_options=jit_options,
               _strict=strict, _inline=staticmethod(InlineOptions(inline)),
               _inline_overloads={}, prefer_literal=prefer_literal,
               metadata=kwargs)
    return type(base)(name, (base,), dct)


class _TemplateTargetHelperMixin(object):
    """Mixin for helper methods that assist with target/registry resolution"""

    def _get_target_registry(self, reason):
        """Returns the registry for the current target.

        Parameters
        ----------
        reason: str
            Reason for the resolution. Expects a noun.
        Returns
        -------
        reg : a registry suitable for the current target.
        """
        from numba.core.target_extension import (_get_local_target_checked,
                                                 dispatcher_registry)
        hwstr = self.metadata.get('target', 'generic')
        target_hw = _get_local_target_checked(self.context, hwstr, reason)
        # Get registry for the current hardware
        disp = dispatcher_registry[target_hw]
        tgtctx = disp.targetdescr.target_context
        # This is all workarounds...
        # The issue is that whilst targets shouldn't care about which registry
        # in which to register lowering implementations, the CUDA target
        # "borrows" implementations from the CPU from specific registries. This
        # means that if some impl is defined via @intrinsic, e.g. numba.*unsafe
        # modules, _AND_ CUDA also makes use of the same impl, then it's
        # required that the registry in use is one that CUDA borrows from. This
        # leads to the following expression where by the CPU builtin_registry is
        # used if it is in the target context as a known registry (i.e. the
        # target installed it) and if it is not then it is assumed that the
        # registries for the target are unbound to any other target and so it's
        # fine to use any of them as a place to put lowering impls.
        #
        # NOTE: This will need subsequently fixing again when targets use solely
        # the extension APIs to describe their implementation. The issue will be
        # that the builtin_registry should contain _just_ the stack allocated
        # implementations and low level target invariant things and should not
        # be modified further. It should be acceptable to remove the `then`
        # branch and just keep the `else`.

        # In case the target has swapped, e.g. cuda borrowing cpu, refresh to
        # populate.
        tgtctx.refresh()
        if builtin_registry in tgtctx._registries:
            reg = builtin_registry
        else:
            # Pick a registry in which to install intrinsics
            registries = iter(tgtctx._registries)
            reg = next(registries)
        return reg


class _IntrinsicTemplate(_TemplateTargetHelperMixin, AbstractTemplate):
    """
    A base class of templates for intrinsic definition
    """

    def generic(self, args, kws):
        """
        Type the intrinsic by the arguments.
        """
        lower_builtin = self._get_target_registry('intrinsic').lower
        cache_key = self.context, args, tuple(kws.items())
        try:
            return self._impl_cache[cache_key]
        except KeyError:
            pass
        result = self._definition_func(self.context, *args, **kws)
        if result is None:
            return
        [sig, imp] = result
        pysig = utils.pysignature(self._definition_func)
        # omit context argument from user function
        parameters = list(pysig.parameters.values())[1:]
        sig = sig.replace(pysig=pysig.replace(parameters=parameters))
        self._impl_cache[cache_key] = sig
        self._overload_cache[sig.args] = imp
        # register the lowering
        lower_builtin(imp, *sig.args)(imp)
        return sig

    def get_impl_key(self, sig):
        """
        Return the key for looking up the implementation for the given
        signature on the target context.
        """
        return self._overload_cache[sig.args]

    def get_template_info(self):
        basepath = os.path.dirname(os.path.dirname(numba.__file__))
        impl = self._definition_func
        code, firstlineno, path = self.get_source_code_info(impl)
        sig = str(utils.pysignature(impl))
        info = {
            'kind': "intrinsic",
            'name': getattr(impl, '__qualname__', impl.__name__),
            'sig': sig,
            'filename': utils.safe_relpath(path, start=basepath),
            'lines': (firstlineno, firstlineno + len(code) - 1),
            'docstring': impl.__doc__
        }
        return info


def make_intrinsic_template(handle, defn, name, *, prefer_literal=False,
                            kwargs=None):
    """
    Make a template class for a intrinsic handle *handle* defined by the
    function *defn*.  The *name* is used for naming the new template class.
    """
    kwargs = MappingProxyType({} if kwargs is None else kwargs)
    base = _IntrinsicTemplate
    name = "_IntrinsicTemplate_%s" % (name)
    dct = dict(key=handle, _definition_func=staticmethod(defn),
               _impl_cache={}, _overload_cache={},
               prefer_literal=prefer_literal, metadata=kwargs)
    return type(base)(name, (base,), dct)


class AttributeTemplate(object):
    def __init__(self, context):
        self.context = context

    def resolve(self, value, attr):
        return self._resolve(value, attr)

    def _resolve(self, value, attr):
        fn = getattr(self, "resolve_%s" % attr, None)
        if fn is None:
            fn = self.generic_resolve
            if fn is NotImplemented:
                if isinstance(value, types.Module):
                    return self.context.resolve_module_constants(value, attr)
                else:
                    return None
            else:
                return fn(value, attr)
        else:
            return fn(value)

    generic_resolve = NotImplemented


class _OverloadAttributeTemplate(_TemplateTargetHelperMixin, AttributeTemplate):
    """
    A base class of templates for @overload_attribute functions.
    """
    is_method = False

    def __init__(self, context):
        super(_OverloadAttributeTemplate, self).__init__(context)
        self.context = context
        self._init_once()

    def _init_once(self):
        cls = type(self)
        attr = cls._attr

        lower_getattr = self._get_target_registry('attribute').lower_getattr

        @lower_getattr(cls.key, attr)
        def getattr_impl(context, builder, typ, value):
            typingctx = context.typing_context
            fnty = cls._get_function_type(typingctx, typ)
            sig = cls._get_signature(typingctx, fnty, (typ,), {})
            call = context.get_function(fnty, sig)
            return call(builder, (value,))

    def _resolve(self, typ, attr):
        if self._attr != attr:
            return None
        fnty = self._get_function_type(self.context, typ)
        sig = self._get_signature(self.context, fnty, (typ,), {})
        # There should only be one template
        for template in fnty.templates:
            self._inline_overloads.update(template._inline_overloads)
        return sig.return_type

    @classmethod
    def _get_signature(cls, typingctx, fnty, args, kws):
        sig = fnty.get_call_type(typingctx, args, kws)
        sig = sig.replace(pysig=utils.pysignature(cls._overload_func))
        return sig

    @classmethod
    def _get_function_type(cls, typingctx, typ):
        return typingctx.resolve_value_type(cls._overload_func)


class _OverloadMethodTemplate(_OverloadAttributeTemplate):
    """
    A base class of templates for @overload_method functions.
    """
    is_method = True

    def _init_once(self):
        """
        Overriding parent definition
        """
        attr = self._attr

        registry = self._get_target_registry('method')

        @registry.lower((self.key, attr), self.key, types.VarArg(types.Any))
        def method_impl(context, builder, sig, args):
            typ = sig.args[0]
            typing_context = context.typing_context
            fnty = self._get_function_type(typing_context, typ)
            sig = self._get_signature(typing_context, fnty, sig.args, {})
            call = context.get_function(fnty, sig)
            # Link dependent library
            context.add_linking_libs(getattr(call, 'libs', ()))
            return call(builder, args)

    def _resolve(self, typ, attr):
        if self._attr != attr:
            return None

        if isinstance(typ, types.TypeRef):
            assert typ == self.key
        elif isinstance(typ, types.Callable):
            assert typ == self.key
        else:
            assert isinstance(typ, self.key)

        class MethodTemplate(AbstractTemplate):
            key = (self.key, attr)
            _inline = self._inline
            _overload_func = staticmethod(self._overload_func)
            _inline_overloads = self._inline_overloads
            prefer_literal = self.prefer_literal

            def generic(_, args, kws):
                args = (typ,) + tuple(args)
                fnty = self._get_function_type(self.context, typ)
                sig = self._get_signature(self.context, fnty, args, kws)
                sig = sig.replace(pysig=utils.pysignature(self._overload_func))
                for template in fnty.templates:
                    self._inline_overloads.update(template._inline_overloads)
                if sig is not None:
                    return sig.as_method()

            def get_template_info(self):
                basepath = os.path.dirname(os.path.dirname(numba.__file__))
                impl = self._overload_func
                code, firstlineno, path = self.get_source_code_info(impl)
                sig = str(utils.pysignature(impl))
                info = {
                    'kind': "overload_method",
                    'name': getattr(impl, '__qualname__', impl.__name__),
                    'sig': sig,
                    'filename': utils.safe_relpath(path, start=basepath),
                    'lines': (firstlineno, firstlineno + len(code) - 1),
                    'docstring': impl.__doc__
                }

                return info

        return types.BoundFunction(MethodTemplate, typ)


def make_overload_attribute_template(typ, attr, overload_func, inline='never',
                                     prefer_literal=False,
                                     base=_OverloadAttributeTemplate,
                                     **kwargs):
    """
    Make a template class for attribute *attr* of *typ* overloaded by
    *overload_func*.
    """
    assert isinstance(typ, types.Type) or issubclass(typ, types.Type)
    name = "OverloadAttributeTemplate_%s_%s" % (typ, attr)
    # Note the implementation cache is subclass-specific
    dct = dict(key=typ, _attr=attr, _impl_cache={},
               _inline=staticmethod(InlineOptions(inline)),
               _inline_overloads={},
               _overload_func=staticmethod(overload_func),
               prefer_literal=prefer_literal,
               metadata=kwargs,
               )
    obj = type(base)(name, (base,), dct)
    return obj


def make_overload_method_template(typ, attr, overload_func, inline,
                                  prefer_literal=False, **kwargs):
    """
    Make a template class for method *attr* of *typ* overloaded by
    *overload_func*.
    """
    return make_overload_attribute_template(
        typ, attr, overload_func, inline=inline,
        base=_OverloadMethodTemplate, prefer_literal=prefer_literal,
        **kwargs,
    )


def bound_function(template_key):
    """
    Wrap an AttributeTemplate resolve_* method to allow it to
    resolve an instance method's signature rather than a instance attribute.
    The wrapped method must return the resolved method's signature
    according to the given self type, args, and keywords.

    It is used thusly:

        class ComplexAttributes(AttributeTemplate):
            @bound_function("complex.conjugate")
            def resolve_conjugate(self, ty, args, kwds):
                return ty

    *template_key* (e.g. "complex.conjugate" above) will be used by the
    target to look up the method's implementation, as a regular function.
    """
    def wrapper(method_resolver):
        @functools.wraps(method_resolver)
        def attribute_resolver(self, ty):
            class MethodTemplate(AbstractTemplate):
                key = template_key

                def generic(_, args, kws):
                    sig = method_resolver(self, ty, args, kws)
                    if sig is not None and sig.recvr is None:
                        sig = sig.replace(recvr=ty)
                    return sig

            return types.BoundFunction(MethodTemplate, ty)
        return attribute_resolver
    return wrapper


# -----------------------------

class Registry(object):
    """
    A registry of typing declarations.  The registry stores such declarations
    for functions, attributes and globals.
    """

    def __init__(self):
        self.functions = []
        self.attributes = []
        self.globals = []

    def register(self, item):
        assert issubclass(item, FunctionTemplate)
        self.functions.append(item)
        return item

    def register_attr(self, item):
        assert issubclass(item, AttributeTemplate)
        self.attributes.append(item)
        return item

    def register_global(self, val=None, typ=None, **kwargs):
        """
        Register the typing of a global value.
        Functional usage with a Numba type::
            register_global(value, typ)

        Decorator usage with a template class::
            @register_global(value, typing_key=None)
            class Template:
                ...
        """
        if typ is not None:
            # register_global(val, typ)
            assert val is not None
            assert not kwargs
            self.globals.append((val, typ))
        else:
            def decorate(cls, typing_key):
                class Template(cls):
                    key = typing_key
                if callable(val):
                    typ = types.Function(Template)
                else:
                    raise TypeError("cannot infer type for global value %r")
                self.globals.append((val, typ))
                return cls

            # register_global(val, typing_key=None)(<template class>)
            assert val is not None
            typing_key = kwargs.pop('typing_key', val)
            assert not kwargs
            if typing_key is val:
                # Check the value is globally reachable, as it is going
                # to be used as the key.
                mod = sys.modules[val.__module__]
                if getattr(mod, val.__name__) is not val:
                    raise ValueError("%r is not globally reachable as '%s.%s'"
                                     % (mod, val.__module__, val.__name__))

            def decorator(cls):
                return decorate(cls, typing_key)
            return decorator


class BaseRegistryLoader(object):
    """
    An incremental loader for a registry.  Each new call to
    new_registrations() will iterate over the not yet seen registrations.

    The reason for this object is multiple:
    - there can be several contexts
    - each context wants to install all registrations
    - registrations can be added after the first installation, so contexts
      must be able to get the "new" installations

    Therefore each context maintains its own loaders for each existing
    registry, without duplicating the registries themselves.
    """

    def __init__(self, registry):
        self._registrations = dict(
            (name, utils.stream_list(getattr(registry, name)))
            for name in self.registry_items)

    def new_registrations(self, name):
        for item in next(self._registrations[name]):
            yield item


class RegistryLoader(BaseRegistryLoader):
    """
    An incremental loader for a typing registry.
    """
    registry_items = ('functions', 'attributes', 'globals')


builtin_registry = Registry()
infer = builtin_registry.register
infer_getattr = builtin_registry.register_attr
infer_global = builtin_registry.register_global
