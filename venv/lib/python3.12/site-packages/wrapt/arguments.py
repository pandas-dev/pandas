"""The inspect.formatargspec() function was dropped in Python 3.11 but we need
it for when constructing signature changing decorators based on result of
inspect.getfullargspec(). The code here implements inspect.formatargspec() based
on Parameter and Signature from inspect module, which were added in Python 3.6.
Thanks to Cyril Jouve for the implementation.
"""

from typing import Any, Callable, List, Mapping, Optional, Sequence, Tuple

try:
    from inspect import Parameter, Signature
except ImportError:
    from inspect import formatargspec  # type: ignore[attr-defined]
else:

    def formatargspec(
        args: List[str],
        varargs: Optional[str] = None,
        varkw: Optional[str] = None,
        defaults: Optional[Tuple[Any, ...]] = None,
        kwonlyargs: Optional[Sequence[str]] = None,
        kwonlydefaults: Optional[Mapping[str, Any]] = None,
        annotations: Mapping[str, Any] = {},
        formatarg: Callable[[str], str] = str,
        formatvarargs: Callable[[str], str] = lambda name: "*" + name,
        formatvarkw: Callable[[str], str] = lambda name: "**" + name,
        formatvalue: Callable[[Any], str] = lambda value: "=" + repr(value),
        formatreturns: Callable[[Any], str] = lambda text: " -> " + text,
        formatannotation: Callable[[Any], str] = lambda annot: " -> " + repr(annot),
    ) -> str:
        if kwonlyargs is None:
            kwonlyargs = ()
        if kwonlydefaults is None:
            kwonlydefaults = {}
        ndefaults = len(defaults) if defaults else 0
        parameters = [
            Parameter(
                arg,
                Parameter.POSITIONAL_OR_KEYWORD,
                default=defaults[i] if defaults and i >= 0 else Parameter.empty,
                annotation=annotations.get(arg, Parameter.empty),
            )
            for i, arg in enumerate(args, ndefaults - len(args))
        ]
        if varargs:
            parameters.append(Parameter(varargs, Parameter.VAR_POSITIONAL))
        parameters.extend(
            Parameter(
                kwonlyarg,
                Parameter.KEYWORD_ONLY,
                default=kwonlydefaults.get(kwonlyarg, Parameter.empty),
                annotation=annotations.get(kwonlyarg, Parameter.empty),
            )
            for kwonlyarg in kwonlyargs
        )
        if varkw:
            parameters.append(Parameter(varkw, Parameter.VAR_KEYWORD))
        return_annotation = annotations.get("return", Signature.empty)
        return str(Signature(parameters, return_annotation=return_annotation))
