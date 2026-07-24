import abc
import collections
import inspect
from collections.abc import Callable, Sequence
from types import CodeType, MethodType
from typing import (
    Any,
    Literal,
    NamedTuple,
    Protocol,
    TypeAlias,
    overload,
    type_check_only,
)
from typing_extensions import (
    ClassVar,
    Generic,
    Never,
    ParamSpec,
    Self,
    TypeVar,
    TypedDict,
    Unpack,
)
from _typeshed import Incomplete, SupportsWrite

from . import base, codegen, compiler, types
from .annotations import pretty_annotate
from .typing import context, templates
from numba import _dispatcher

_T = TypeVar("_T")
_T_co = TypeVar("_T_co", default=Any, covariant=True)
_FuncT_co = TypeVar(
    "_FuncT_co",
    bound=Callable[..., object],
    default=Callable[..., Any],
    covariant=True,
)
_ParamsT = ParamSpec("_ParamsT")

_InputSignature: TypeAlias = tuple[types.Type, ...]
_ToSignature: TypeAlias = str | _InputSignature | templates.Signature

_WrapperKind: TypeAlias = Literal["python", "cfunc"]

@type_check_only
class _CFGKwargs(TypedDict, total=False):
    filename: str
    view: bool
    highlight: bool | set[str] | dict[str, bool]
    interleave: bool | set[str] | dict[str, bool]
    strip_ir: bool
    show_key: bool
    fontsize: int

# represents `numba.misc.inspection.disassemble_elf_to_cfg.<locals>.DisasmCFG`
@type_check_only
class _DisasmCFG(Protocol):
    def _repr_svg_(self) -> str: ...
    def __repr__(self) -> str: ...

###

# undocumented
class OmittedArg(Generic[_T_co]):
    value: _T_co

    def __init__(self, value: _T_co) -> None: ...
    @property
    def _numba_type_(self) -> types.Omitted: ...

# undocumented
class CompilingCounter:
    counter: int

    def __init__(self) -> None: ...
    def __enter__(self) -> None: ...
    def __exit__(self, *args: Any, **kwargs: Any) -> None: ...
    def __bool__(self) -> bool: ...

# undocumented
class _CompileStats(NamedTuple):
    cache_path: str | None
    cache_hits: collections.Counter[str]
    cache_misses: collections.Counter[str]

# NOTE: `_DispatcherBase` is not a generic type at runtime.
# undocumented
class _DispatcherBase(_dispatcher.Dispatcher[_FuncT_co], Generic[_FuncT_co]):
    _fold_args: ClassVar[bool]  # abstract
    __numba__: ClassVar[Literal["py_func"]] = "py_func"

    overloads: collections.OrderedDict[tuple[types.Type, ...], compiler.CompileResult]
    py_func: _FuncT_co

    __code__: CodeType
    __defaults__: tuple[Any, ...]
    doc: str | None

    _compiling_counter: CompilingCounter
    _enable_sysmon: bool

    def __init__(
        self,
        arg_count: int,
        py_func: _FuncT_co,
        pysig: inspect.Signature,
        can_fallback: bool,
        exact_match_required: bool,
    ) -> None: ...

    #
    @property
    def signatures(self) -> list[_InputSignature]: ...
    @property
    def nopython_signatures(self) -> list[templates.Signature]: ...

    #
    def disable_compile(self, val: bool = True) -> None: ...
    def add_overload(self, cres: compiler.CompileResult) -> None: ...
    def fold_argument_types(
        self, args: Sequence[_T], kws: Sequence[tuple[str, _T]] | dict[str, _T]
    ) -> tuple[inspect.Signature, tuple[_T, ...]]: ...
    def get_call_template(
        self, args: Sequence[_T], kws: Sequence[tuple[str, _T]] | dict[str, _T]
    ) -> tuple[
        type[templates.CallableTemplate],
        inspect.Signature,
        tuple[_T, ...],
        dict[Never, Never],
    ]: ...
    def get_overload(
        self: _DispatcherBase[Callable[_ParamsT, _T]], sig: _ToSignature
    ) -> Callable[_ParamsT, _T]: ...

    #
    @property
    def is_compiling(self) -> CompilingCounter: ...

    #
    @overload  # signature=None (default)
    def inspect_llvm(self, signature: None = None) -> dict[_InputSignature, str]: ...
    @overload  # signature=<given>
    def inspect_llvm(self, signature: _InputSignature) -> str: ...

    #
    @overload  # signature=None (default)
    def inspect_asm(self, signature: None = None) -> dict[_InputSignature, str]: ...
    @overload  # signature=<given>
    def inspect_asm(self, signature: _InputSignature) -> str: ...

    #
    @overload  # pretty=False (default)
    def inspect_types(
        self,
        file: SupportsWrite[str] | None = None,
        signature: _InputSignature | None = None,
        pretty: Literal[False] = False,
        style: str = "default",
        **kwargs: Never,
    ) -> None: ...
    @overload  # pretty=True (positional)
    def inspect_types(
        self,
        file: None,
        signature: _InputSignature | None,
        pretty: Literal[True],
        style: str = "default",
        **kwargs: Never,
    ) -> pretty_annotate.Annotate: ...
    @overload  # pretty=True (keyword)
    def inspect_types(
        self,
        file: None = None,
        signature: _InputSignature | None = None,
        *,
        pretty: Literal[True],
        style: str = "default",
    ) -> pretty_annotate.Annotate: ...

    #
    @overload  # signature=None (default)
    def inspect_cfg(
        self,
        signature: None = None,
        show_wrapper: _WrapperKind | None = None,
        **kwargs: Unpack[_CFGKwargs],
    ) -> dict[_InputSignature, codegen._CFG]: ...
    @overload  # signature=<given>
    def inspect_cfg(
        self,
        signature: _InputSignature,
        show_wrapper: _WrapperKind | None = None,
        **kwargs: Unpack[_CFGKwargs],
    ) -> codegen._CFG: ...

    #
    @overload  # signature=None (default)
    def inspect_disasm_cfg(
        self, signature: None = None
    ) -> dict[_InputSignature, _DisasmCFG]: ...
    @overload  # signature=<given>
    def inspect_disasm_cfg(self, signature: _InputSignature) -> _DisasmCFG: ...

    #
    def get_annotation_info(
        self, signature: _InputSignature | None = None
    ) -> collections.OrderedDict[tuple[str, str], dict[str, Incomplete]]: ...

    #
    def typeof_pyval(self, val: object) -> types.Type: ...

class Dispatcher(_DispatcherBase[_FuncT_co], Generic[_FuncT_co], metaclass=abc.ABCMeta):
    _fold_args: ClassVar[bool] = True
    __numba__: ClassVar[Literal["py_func"]] = "py_func"

    typingctx: context.Context
    targetctx: base.BaseContext
    targetoptions: dict[str, Any]
    locals: dict[str, Any]

    def __init__(
        self,
        py_func: _FuncT_co,
        locals: dict[str, Any] | None = None,
        targetoptions: dict[str, Any] | None = None,
        pipeline_class: type[compiler.CompilerBase] = compiler.Compiler,
    ) -> None: ...
    def __reduce__(self) -> tuple[Any, Any]: ...
    def dump(self, tab: str = "") -> None: ...

    #
    @property
    def _numba_type_(self) -> types.Dispatcher: ...

    #
    def enable_caching(self) -> None: ...

    #
    @overload  # obj=None
    def __get__(self, obj: None, objtype: type | None = None) -> Self: ...
    @overload  # obj=<given>
    def __get__(self, obj: object, objtype: type | None = None) -> MethodType: ...

    #
    def compile(
        self: Dispatcher[Callable[_ParamsT, _T]], sig: _ToSignature
    ) -> Callable[_ParamsT, _T]: ...
    def get_compile_result(
        self, sig: templates.Signature
    ) -> compiler.CompileResult: ...
    def recompile(self) -> None: ...

    #
    @property
    def stats(self) -> _CompileStats: ...

    #
    def parallel_diagnostics(
        self, signature: _InputSignature | None = None, level: int = 1
    ) -> None: ...

    #
    @overload  # signature=None (default)
    def get_metadata(
        self, signature: None = None
    ) -> dict[_InputSignature, dict[str, Incomplete]]: ...
    @overload  # signature=<given>
    def get_metadata(self, signature: _InputSignature) -> dict[str, Incomplete]: ...

    #
    def get_function_type(self) -> types.FunctionType | None: ...


# TODO
LiftedCode: Incomplete
LiftedLoop: Incomplete
LiftedWith: Incomplete
ObjModeLiftedWith: Incomplete
