"""Tools for inspecting Python objects.

Uses syntax highlighting for presenting the various information elements.

Similar in spirit to the inspect module, but all calls take a name argument to
reference the name under which an object is being read.
"""

# Copyright (c) IPython Development Team.
# Distributed under the terms of the Modified BSD License.

__all__ = ["Inspector"]

# stdlib modules
from dataclasses import dataclass
from inspect import signature
from textwrap import dedent
import ast
import html
import inspect
import io as stdlib_io
import linecache
import os
import types
import warnings
from pygments.token import Token


from typing import (
    cast,
    Any,
    Optional,
    Dict,
    Union,
    List,
    TypedDict,
    TypeAlias,
    Tuple,
)

import traitlets
from traitlets.config import Configurable

# IPython's own
from IPython.core import page
from IPython.lib.pretty import pretty
from IPython.testing.skipdoctest import skip_doctest
from IPython.utils import PyColorize, openpy
from IPython.utils.dir2 import safe_hasattr
from IPython.utils.path import compress_user
from IPython.utils.text import indent
from IPython.utils.wildcard import list_namespace, typestr2type
from IPython.utils.decorators import undoc

from pygments import highlight
from pygments.lexers import PythonLexer
from pygments.formatters import HtmlFormatter

HOOK_NAME = "__custom_documentations__"


UnformattedBundle: TypeAlias = Dict[str, List[Tuple[str, str]]]  # List of (title, body)
Bundle: TypeAlias = Dict[str, str]


@dataclass
class OInfo:
    ismagic: bool
    isalias: bool
    found: bool
    namespace: Optional[str]
    parent: Any
    obj: Any

    def get(self, field):
        """Get a field from the object for backward compatibility with before 8.12

        see https://github.com/h5py/h5py/issues/2253
        """
        # We need to deprecate this at some point, but the warning will show in completion.
        # Let's comment this for now and uncomment end of 2023 ish
        # Jan 2025: decomenting for IPython 9.0
        warnings.warn(
            f"OInfo dataclass with fields access since IPython 8.12 please use OInfo.{field} instead."
            "OInfo used to be a dict but a dataclass provide static fields verification with mypy."
            "This warning and backward compatibility `get()` method were added in 8.13.",
            DeprecationWarning,
            stacklevel=2,
        )
        return getattr(self, field)


def pylight(code):
    return highlight(code, PythonLexer(), HtmlFormatter(noclasses=True))

# builtin docstrings to ignore
_func_call_docstring = types.FunctionType.__call__.__doc__
_object_init_docstring = object.__init__.__doc__
_builtin_type_docstrings = {
    inspect.getdoc(t) for t in (types.ModuleType, types.MethodType,
                                types.FunctionType, property)
}

_builtin_func_type = type(all)
_builtin_meth_type = type(str.upper)  # Bound methods have the same type as builtin functions
#****************************************************************************
# Builtin color schemes


#****************************************************************************
# Auxiliary functions and objects


class InfoDict(TypedDict):
    type_name: Optional[str]
    base_class: Optional[str]
    string_form: Optional[str]
    namespace: Optional[str]
    length: Optional[str]
    file: Optional[str]
    definition: Optional[str]
    docstring: Optional[str]
    source: Optional[str]
    init_definition: Optional[str]
    class_docstring: Optional[str]
    init_docstring: Optional[str]
    call_def: Optional[str]
    call_docstring: Optional[str]
    subclasses: Optional[str]
    # These won't be printed but will be used to determine how to
    # format the object
    ismagic: bool
    isalias: bool
    isclass: bool
    found: bool
    name: str


_info_fields = list(InfoDict.__annotations__.keys())


def __getattr__(name):
    if name == "info_fields":
        warnings.warn(
            "IPython.core.oinspect's `info_fields` is considered for deprecation and may be removed in the Future. ",
            DeprecationWarning,
            stacklevel=2,
        )
        return _info_fields

    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")


@dataclass
class InspectorHookData:
    """Data passed to the mime hook"""

    obj: Any
    info: Optional[OInfo]
    info_dict: InfoDict
    detail_level: int
    omit_sections: list[str]


@undoc
def object_info(
    *,
    name: str,
    found: bool,
    isclass: bool = False,
    isalias: bool = False,
    ismagic: bool = False,
    **kw,
) -> InfoDict:
    """Make an object info dict with all fields present."""
    infodict = dict(kw)
    infodict.update({k: None for k in _info_fields if k not in infodict})
    infodict["name"] = name  # type: ignore
    infodict["found"] = found  # type: ignore
    infodict["isclass"] = isclass  # type: ignore
    infodict["isalias"] = isalias  # type: ignore
    infodict["ismagic"] = ismagic  # type: ignore

    return InfoDict(**infodict)  # type:ignore


def get_encoding(obj):
    """Get encoding for python source file defining obj

    Returns None if obj is not defined in a sourcefile.
    """
    ofile = find_file(obj)
    # run contents of file through pager starting at line where the object
    # is defined, as long as the file isn't binary and is actually on the
    # filesystem.
    if ofile is None:
        return None
    elif ofile.endswith(('.so', '.dll', '.pyd')):
        return None
    elif not os.path.isfile(ofile):
        return None
    else:
        # Print only text files, not extension binaries.  Note that
        # getsourcelines returns lineno with 1-offset and page() uses
        # 0-offset, so we must adjust.
        with stdlib_io.open(ofile, 'rb') as buffer:   # Tweaked to use io.open for Python 2
            encoding, _lines = openpy.detect_encoding(buffer.readline)
        return encoding


def getdoc(obj) -> Union[str, None]:
    """Stable wrapper around inspect.getdoc.

    This can't crash because of attribute problems.

    It also attempts to call a getdoc() method on the given object.  This
    allows objects which provide their docstrings via non-standard mechanisms
    (like Pyro proxies) to still be inspected by ipython's ? system.
    """
    # Allow objects to offer customized documentation via a getdoc method:
    try:
        ds = obj.getdoc()
    except Exception:
        pass
    else:
        if isinstance(ds, str):
            return inspect.cleandoc(ds)
    docstr = inspect.getdoc(obj)
    return docstr


def getsource(obj, oname='') -> Union[str,None]:
    """Wrapper around inspect.getsource.

    This can be modified by other projects to provide customized source
    extraction.

    Parameters
    ----------
    obj : object
        an object whose source code we will attempt to extract
    oname : str
        (optional) a name under which the object is known

    Returns
    -------
    src : unicode or None

    """

    if isinstance(obj, property):
        sources = []
        for attrname in ['fget', 'fset', 'fdel']:
            fn = getattr(obj, attrname)
            if fn is not None:
                oname_prefix = ('%s.' % oname) if oname else ''
                sources.append(''.join(('# ', oname_prefix, attrname)))
                if inspect.isfunction(fn):
                    _src = getsource(fn)
                    if _src:
                        # assert _src is not None, "please mypy"
                        sources.append(dedent(_src))
                else:
                    # Default str/repr only prints function name,
                    # pretty.pretty prints module name too.
                    sources.append(
                        '%s%s = %s\n' % (oname_prefix, attrname, pretty(fn))
                    )
        if sources:
            return '\n'.join(sources)
        else:
            return None

    else:
        # Get source for non-property objects.

        obj = _get_wrapped(obj)

        try:
            src = inspect.getsource(obj)
        except TypeError:
            # The object itself provided no meaningful source, try looking for
            # its class definition instead.
            try:
                src = inspect.getsource(obj.__class__)
            except (OSError, TypeError):
                return None
        except OSError:
            return None

        return src


def is_simple_callable(obj):
    """True if obj is a function ()"""
    return (inspect.isfunction(obj) or inspect.ismethod(obj) or \
            isinstance(obj, _builtin_func_type) or isinstance(obj, _builtin_meth_type))

def _get_wrapped(obj):
    """Get the original object if wrapped in one or more @decorators

    Some objects automatically construct similar objects on any unrecognised
    attribute access (e.g. unittest.mock.call). To protect against infinite loops,
    this will arbitrarily cut off after 100 levels of obj.__wrapped__
    attribute access. --TK, Jan 2016
    """
    orig_obj = obj
    i = 0
    while safe_hasattr(obj, '__wrapped__'):
        obj = obj.__wrapped__
        i += 1
        if i > 100:
            # __wrapped__ is probably a lie, so return the thing we started with
            return orig_obj
    return obj

def find_file(obj) -> Optional[str]:
    """Find the absolute path to the file where an object was defined.

    This is essentially a robust wrapper around `inspect.getabsfile`.

    Returns None if no file can be found.

    Parameters
    ----------
    obj : any Python object

    Returns
    -------
    fname : str
        The absolute path to the file where the object was defined.
    """
    obj = _get_wrapped(obj)

    fname: Optional[str] = None
    try:
        fname = inspect.getabsfile(obj)
    except TypeError:
        # For an instance, the file that matters is where its class was
        # declared.
        try:
            fname = inspect.getabsfile(obj.__class__)
        except (OSError, TypeError):
            # Can happen for builtins
            pass
    except OSError:
        pass

    return fname


def find_source_lines(obj):
    """Find the line number in a file where an object was defined.

    This is essentially a robust wrapper around `inspect.getsourcelines`.

    Returns None if no file can be found.

    Parameters
    ----------
    obj : any Python object

    Returns
    -------
    lineno : int
        The line number where the object definition starts.
    """
    obj = _get_wrapped(obj)

    try:
        lineno = inspect.getsourcelines(obj)[1]
    except TypeError:
        # For instances, try the class object like getsource() does
        try:
            lineno = inspect.getsourcelines(obj.__class__)[1]
        except (OSError, TypeError):
            return None
    except OSError:
        return None

    return lineno


_sentinel = object()


class Inspector(Configurable):
    mime_hooks = traitlets.Dict(
        config=True,
        help="dictionary of mime to callable to add information into help mimebundle dict",
    ).tag(config=True)

    _theme_name: str

    def __init__(
        self,
        *,
        theme_name: str,
        str_detail_level=0,
        parent=None,
        config=None,
    ):
        if theme_name in ["Linux", "LightBG", "Neutral", "NoColor"]:
            warnings.warn(
                f"Theme names and color schemes are lowercase in IPython 9.0 use {theme_name.lower()} instead",
                DeprecationWarning,
                stacklevel=2,
            )
            theme_name = theme_name.lower()
        self._theme_name = theme_name
        super(Inspector, self).__init__(parent=parent, config=config)
        self.parser = PyColorize.Parser(out="str", theme_name=theme_name)
        self.str_detail_level = str_detail_level
        self.set_theme_name(theme_name)

    def format(self, *args, **kwargs):
        return self.parser.format(*args, **kwargs)

    def _getdef(self,obj,oname='') -> Union[str,None]:
        """Return the call signature for any callable object.

        If any exception is generated, None is returned instead and the
        exception is suppressed."""
        if not callable(obj):
            return None
        try:
            return _render_signature(signature(obj), oname)
        except:
            return None

    def __head(self, h: str) -> str:
        """Return a header string with proper colors."""
        return PyColorize.theme_table[self._theme_name].format([(Token.Header, h)])

    def set_theme_name(self, name: str):
        assert name == name.lower()
        assert name in PyColorize.theme_table.keys()
        self._theme_name = name
        self.parser.theme_name = name

    def set_active_scheme(self, scheme: str):
        warnings.warn(
            "set_active_scheme is deprecated and replaced by set_theme_name as of IPython 9.0",
            DeprecationWarning,
            stacklevel=2,
        )
        assert scheme == scheme.lower()
        if scheme is not None and self._theme_name != scheme:
            self._theme_name = scheme
            self.parser.theme_name = scheme

    def noinfo(self, msg, oname):
        """Generic message when no information is found."""
        print('No %s found' % msg, end=' ')
        if oname:
            print('for %s' % oname)
        else:
            print()

    def pdef(self, obj, oname=''):
        """Print the call signature for any callable object.

        If the object is a class, print the constructor information."""

        if not callable(obj):
            print('Object is not callable.')
            return

        header = ''

        if inspect.isclass(obj):
            header = self.__head('Class constructor information:\n')


        output = self._getdef(obj,oname)
        if output is None:
            self.noinfo('definition header',oname)
        else:
            print(header,self.format(output), end=' ')

    # In Python 3, all classes are new-style, so they all have __init__.
    @skip_doctest
    def pdoc(self, obj, oname='', formatter=None):
        """Print the docstring for any object.

        Optional:
        -formatter: a function to run the docstring through for specially
        formatted docstrings.

        Examples
        --------
        In [1]: class NoInit:
           ...:     pass

        In [2]: class NoDoc:
           ...:     def __init__(self):
           ...:         pass

        In [3]: %pdoc NoDoc
        No documentation found for NoDoc

        In [4]: %pdoc NoInit
        No documentation found for NoInit

        In [5]: obj = NoInit()

        In [6]: %pdoc obj
        No documentation found for obj

        In [5]: obj2 = NoDoc()

        In [6]: %pdoc obj2
        No documentation found for obj2
        """

        lines = []
        ds = getdoc(obj)
        if formatter:
            ds = formatter(ds).get('plain/text', ds)
        if ds:
            lines.append(self.__head("Class docstring:"))
            lines.append(indent(ds))
        if inspect.isclass(obj) and hasattr(obj, '__init__'):
            init_ds = getdoc(obj.__init__)
            if init_ds is not None:
                lines.append(self.__head("Init docstring:"))
                lines.append(indent(init_ds))
        elif hasattr(obj,'__call__'):
            call_ds = getdoc(obj.__call__)
            if call_ds:
                lines.append(self.__head("Call docstring:"))
                lines.append(indent(call_ds))

        if not lines:
            self.noinfo('documentation',oname)
        else:
            page.page('\n'.join(lines))

    def psource(self, obj, oname=''):
        """Print the source code for an object."""

        # Flush the source cache because inspect can return out-of-date source
        linecache.checkcache()
        try:
            src = getsource(obj, oname=oname)
        except Exception:
            src = None

        if src is None:
            self.noinfo('source', oname)
        else:
            page.page(self.format(src))

    def pfile(self, obj, oname=''):
        """Show the whole file where an object was defined."""

        lineno = find_source_lines(obj)
        if lineno is None:
            self.noinfo('file', oname)
            return

        ofile = find_file(obj)
        # run contents of file through pager starting at line where the object
        # is defined, as long as the file isn't binary and is actually on the
        # filesystem.
        if ofile is None:
            print("Could not find file for object")
        elif ofile.endswith((".so", ".dll", ".pyd")):
            print("File %r is binary, not printing." % ofile)
        elif not os.path.isfile(ofile):
            print('File %r does not exist, not printing.' % ofile)
        else:
            # Print only text files, not extension binaries.  Note that
            # getsourcelines returns lineno with 1-offset and page() uses
            # 0-offset, so we must adjust.
            page.page(self.format(openpy.read_py_file(ofile, skip_encoding_cookie=False)), lineno - 1)


    def _mime_format(self, text:str, formatter=None) -> dict:
        """Return a mime bundle representation of the input text.

        - if `formatter` is None, the returned mime bundle has
           a ``text/plain`` field, with the input text.
           a ``text/html`` field with a ``<pre>`` tag containing the input text.

        - if ``formatter`` is not None, it must be a callable transforming the
          input text into a mime bundle. Default values for ``text/plain`` and
          ``text/html`` representations are the ones described above.

        Note:

        Formatters returning strings are supported but this behavior is deprecated.

        """
        defaults = {
            "text/plain": text,
            "text/html": f"<pre>{html.escape(text)}</pre>",
        }

        if formatter is None:
            return defaults
        else:
            formatted = formatter(text)

            if not isinstance(formatted, dict):
                # Handle the deprecated behavior of a formatter returning
                # a string instead of a mime bundle.
                return {"text/plain": formatted, "text/html": f"<pre>{formatted}</pre>"}

            else:
                return dict(defaults, **formatted)

    def format_mime(self, bundle: UnformattedBundle) -> Bundle:
        """Format a mimebundle being created by _make_info_unformatted into a real mimebundle"""
        # Format text/plain mimetype
        assert isinstance(bundle["text/plain"], list)
        for item in bundle["text/plain"]:
            assert isinstance(item, tuple)

        new_b: Bundle = {}
        lines = []
        _len = max(len(h) for h, _ in bundle["text/plain"])

        for head, body in bundle["text/plain"]:
            body = body.strip("\n")
            delim = "\n" if "\n" in body else " "
            lines.append(
                f"{self.__head(head+':')}{(_len - len(head))*' '}{delim}{body}"
            )

        new_b["text/plain"] = "\n".join(lines)

        if "text/html" in bundle:
            assert isinstance(bundle["text/html"], list)
            for item in bundle["text/html"]:
                assert isinstance(item, tuple)
            # Format the text/html mimetype
            if isinstance(bundle["text/html"], (list, tuple)):
                # bundle['text/html'] is a list of (head, formatted body) pairs
                new_b["text/html"] = "\n".join(
                    f"<h1>{head}</h1>\n{body}" for (head, body) in bundle["text/html"]
                )

        for k in bundle.keys():
            if k in ("text/html", "text/plain"):
                continue
            else:
                new_b[k] = bundle[k]  # type:ignore
        return new_b

    def _append_info_field(
        self,
        bundle: UnformattedBundle,
        title: str,
        key: str,
        info,
        omit_sections: List[str],
        formatter,
    ):
        """Append an info value to the unformatted mimebundle being constructed by _make_info_unformatted"""
        if title in omit_sections or key in omit_sections:
            return
        field = info[key]
        if field is not None:
            formatted_field = self._mime_format(field, formatter)
            bundle["text/plain"].append((title, formatted_field["text/plain"]))
            bundle["text/html"].append((title, formatted_field["text/html"]))

    def _make_info_unformatted(
        self, obj, info, formatter, detail_level, omit_sections
    ) -> UnformattedBundle:
        """Assemble the mimebundle as unformatted lists of information"""
        bundle: UnformattedBundle = {
            "text/plain": [],
            "text/html": [],
        }

        # A convenience function to simplify calls below
        def append_field(
            bundle: UnformattedBundle, title: str, key: str, formatter=None
        ):
            self._append_info_field(
                bundle,
                title=title,
                key=key,
                info=info,
                omit_sections=omit_sections,
                formatter=formatter,
            )

        def code_formatter(text) -> Bundle:
            return {
                'text/plain': self.format(text),
                'text/html': pylight(text)
            }

        if info["isalias"]:
            append_field(bundle, "Repr", "string_form")

        elif info['ismagic']:
            if detail_level > 0:
                append_field(bundle, "Source", "source", code_formatter)
            else:
                append_field(bundle, "Docstring", "docstring", formatter)
            append_field(bundle, "File", "file")

        elif info['isclass'] or is_simple_callable(obj):
            # Functions, methods, classes
            append_field(bundle, "Signature", "definition", code_formatter)
            append_field(bundle, "Init signature", "init_definition", code_formatter)
            append_field(bundle, "Docstring", "docstring", formatter)
            if detail_level > 0 and info["source"]:
                append_field(bundle, "Source", "source", code_formatter)
            else:
                append_field(bundle, "Init docstring", "init_docstring", formatter)

            append_field(bundle, "File", "file")
            append_field(bundle, "Type", "type_name")
            append_field(bundle, "Subclasses", "subclasses")

        else:
            # General Python objects
            append_field(bundle, "Signature", "definition", code_formatter)
            append_field(bundle, "Call signature", "call_def", code_formatter)
            append_field(bundle, "Type", "type_name")
            append_field(bundle, "String form", "string_form")

            # Namespace
            if info["namespace"] != "Interactive":
                append_field(bundle, "Namespace", "namespace")

            append_field(bundle, "Length", "length")
            append_field(bundle, "File", "file")

            # Source or docstring, depending on detail level and whether
            # source found.
            if detail_level > 0 and info["source"]:
                append_field(bundle, "Source", "source", code_formatter)
            else:
                append_field(bundle, "Docstring", "docstring", formatter)

            append_field(bundle, "Class docstring", "class_docstring", formatter)
            append_field(bundle, "Init docstring", "init_docstring", formatter)
            append_field(bundle, "Call docstring", "call_docstring", formatter)
        return bundle


    def _get_info(
        self,
        obj: Any,
        oname: str = "",
        formatter=None,
        info: Optional[OInfo] = None,
        detail_level: int = 0,
        omit_sections: Union[List[str], Tuple[()]] = (),
    ) -> Bundle:
        """Retrieve an info dict and format it.

        Parameters
        ----------
        obj : any
            Object to inspect and return info from
        oname : str (default: ''):
            Name of the variable pointing to `obj`.
        formatter : callable
        info
            already computed information
        detail_level : integer
            Granularity of detail level, if set to 1, give more information.
        omit_sections : list[str]
            Titles or keys to omit from output (can be set, tuple, etc., anything supporting `in`)
        """

        info_dict = self.info(obj, oname=oname, info=info, detail_level=detail_level)
        omit_sections = list(omit_sections)

        bundle = self._make_info_unformatted(
            obj,
            info_dict,
            formatter,
            detail_level=detail_level,
            omit_sections=omit_sections,
        )
        if self.mime_hooks:
            hook_data = InspectorHookData(
                obj=obj,
                info=info,
                info_dict=info_dict,
                detail_level=detail_level,
                omit_sections=omit_sections,
            )
            for key, hook in self.mime_hooks.items():  # type:ignore
                required_parameters = [
                    parameter
                    for parameter in inspect.signature(hook).parameters.values()
                    if parameter.default != inspect.Parameter.default
                ]
                if len(required_parameters) == 1:
                    res = hook(hook_data)
                else:
                    warnings.warn(
                        "MIME hook format changed in IPython 8.22; hooks should now accept"
                        " a single parameter (InspectorHookData); support for hooks requiring"
                        " two-parameters (obj and info) will be removed in a future version",
                        DeprecationWarning,
                        stacklevel=2,
                    )
                    res = hook(obj, info)
                if res is not None:
                    bundle[key] = res
        return self.format_mime(bundle)

    def pinfo(
        self,
        obj,
        oname="",
        formatter=None,
        info: Optional[OInfo] = None,
        detail_level=0,
        enable_html_pager=True,
        omit_sections=(),
    ):
        """Show detailed information about an object.

        Optional arguments:

        - oname: name of the variable pointing to the object.

        - formatter: callable (optional)
              A special formatter for docstrings.

              The formatter is a callable that takes a string as an input
              and returns either a formatted string or a mime type bundle
              in the form of a dictionary.

              Although the support of custom formatter returning a string
              instead of a mime type bundle is deprecated.

        - info: a structure with some information fields which may have been
          precomputed already.

        - detail_level: if set to 1, more information is given.

        - omit_sections: set of section keys and titles to omit
        """
        assert info is not None
        info_b: Bundle = self._get_info(
            obj, oname, formatter, info, detail_level, omit_sections=omit_sections
        )
        if not enable_html_pager:
            del info_b["text/html"]
        page.page(info_b)

    def info(self, obj, oname="", info=None, detail_level=0) -> InfoDict:
        """Compute a dict with detailed information about an object.

        Parameters
        ----------
        obj : any
            An object to find information about
        oname : str (default: '')
            Name of the variable pointing to `obj`.
        info : (default: None)
            A struct (dict like with attr access) with some information fields
            which may have been precomputed already.
        detail_level : int (default:0)
            If set to 1, more information is given.

        Returns
        -------
        An object info dict with known fields from `info_fields` (see `InfoDict`).
        """

        if info is None:
            ismagic = False
            isalias = False
            ospace = ''
        else:
            ismagic = info.ismagic
            isalias = info.isalias
            ospace = info.namespace

        # Get docstring, special-casing aliases:
        att_name = oname.split(".")[-1]
        parents_docs = None
        prelude = ""
        if info and info.parent is not None and hasattr(info.parent, HOOK_NAME):
            parents_docs_dict = getattr(info.parent, HOOK_NAME)
            parents_docs = parents_docs_dict.get(att_name, None)
        out: InfoDict = cast(
            InfoDict,
            {
                **{field: None for field in _info_fields},
                **{
                    "name": oname,
                    "found": True,
                    "isalias": isalias,
                    "ismagic": ismagic,
                    "subclasses": None,
                },
            },
        )

        if parents_docs:
            ds = parents_docs
        elif isalias:
            if not callable(obj):
                try:
                    ds = "Alias to the system command:\n  %s" % obj[1]
                except:
                    ds = "Alias: " + str(obj)
            else:
                ds = "Alias to " + str(obj)
                if obj.__doc__:
                    ds += "\nDocstring:\n" + obj.__doc__
        else:
            ds_or_None = getdoc(obj)
            if ds_or_None is None:
                ds = '<no docstring>'
            else:
                ds = ds_or_None

        ds = prelude + ds

        # store output in a dict, we initialize it here and fill it as we go

        string_max = 200 # max size of strings to show (snipped if longer)
        shalf = int((string_max - 5) / 2)

        if ismagic:
            out['type_name'] = 'Magic function'
        elif isalias:
            out['type_name'] = 'System alias'
        else:
            out['type_name'] = type(obj).__name__

        try:
            bclass = obj.__class__
            out['base_class'] = str(bclass)
        except:
            pass

        # String form, but snip if too long in ? form (full in ??)
        if detail_level >= self.str_detail_level:
            try:
                ostr = str(obj)
                if not detail_level and len(ostr) > string_max:
                    ostr = ostr[:shalf] + ' <...> ' + ostr[-shalf:]
                    # TODO: `'string_form'.expandtabs()` seems wrong, but
                    # it was (nearly) like this since the first commit ever.
                    ostr = ("\n" + " " * len("string_form".expandtabs())).join(
                        q.strip() for q in ostr.split("\n")
                    )
                out["string_form"] = ostr
            except:
                pass

        if ospace:
            out['namespace'] = ospace

        # Length (for strings and lists)
        try:
            out['length'] = str(len(obj))
        except Exception:
            pass

        # Filename where object was defined
        binary_file = False
        fname = find_file(obj)
        if fname is None:
            # if anything goes wrong, we don't want to show source, so it's as
            # if the file was binary
            binary_file = True
        else:
            if fname.endswith(('.so', '.dll', '.pyd')):
                binary_file = True
            elif fname.endswith('<string>'):
                fname = 'Dynamically generated function. No source code available.'
            out['file'] = compress_user(fname)

        # Original source code for a callable, class or property.
        if detail_level:
            # Flush the source cache because inspect can return out-of-date
            # source
            linecache.checkcache()
            try:
                if isinstance(obj, property) or not binary_file:
                    src = getsource(obj, oname)
                    if src is not None:
                        src = src.rstrip()
                    out['source'] = src

            except Exception:
                pass

        # Add docstring only if no source is to be shown (avoid repetitions).
        if ds and not self._source_contains_docstring(out.get('source'), ds):
            out['docstring'] = ds

        # Constructor docstring for classes
        if inspect.isclass(obj):
            out['isclass'] = True

            # get the init signature:
            try:
                init_def = self._getdef(obj, oname)
            except AttributeError:
                init_def = None

            # get the __init__ docstring
            try:
                obj_init = obj.__init__
            except AttributeError:
                init_ds = None
            else:
                if init_def is None:
                    # Get signature from init if top-level sig failed.
                    # Can happen for built-in types (list, etc.).
                    try:
                        init_def = self._getdef(obj_init, oname)
                    except AttributeError:
                        pass
                init_ds = getdoc(obj_init)
                # Skip Python's auto-generated docstrings
                if init_ds == _object_init_docstring:
                    init_ds = None

            if init_def:
                out['init_definition'] = init_def

            if init_ds:
                out['init_docstring'] = init_ds

            names = [sub.__name__ for sub in type.__subclasses__(obj)]
            if len(names) < 10:
                all_names = ', '.join(names)
            else:
                all_names = ', '.join(names[:10]+['...'])
            out['subclasses'] = all_names
        # and class docstring for instances:
        else:
            # reconstruct the function definition and print it:
            defln = self._getdef(obj, oname)
            if defln:
                out['definition'] = defln

            # First, check whether the instance docstring is identical to the
            # class one, and print it separately if they don't coincide.  In
            # most cases they will, but it's nice to print all the info for
            # objects which use instance-customized docstrings.
            if ds:
                try:
                    cls = getattr(obj,'__class__')
                except:
                    class_ds = None
                else:
                    class_ds = getdoc(cls)
                # Skip Python's auto-generated docstrings
                if class_ds in _builtin_type_docstrings:
                    class_ds = None
                if class_ds and ds != class_ds:
                    out['class_docstring'] = class_ds

            # Next, try to show constructor docstrings
            try:
                init_ds = getdoc(obj.__init__)
                # Skip Python's auto-generated docstrings
                if init_ds == _object_init_docstring:
                    init_ds = None
            except AttributeError:
                init_ds = None
            if init_ds:
                out['init_docstring'] = init_ds

            # Call form docstring for callable instances
            if safe_hasattr(obj, '__call__') and not is_simple_callable(obj):
                call_def = self._getdef(obj.__call__, oname)
                if call_def and (call_def != out.get('definition')):
                    # it may never be the case that call def and definition differ,
                    # but don't include the same signature twice
                    out['call_def'] = call_def
                call_ds = getdoc(obj.__call__)
                # Skip Python's auto-generated docstrings
                if call_ds == _func_call_docstring:
                    call_ds = None
                if call_ds:
                    out['call_docstring'] = call_ds

        return out

    @staticmethod
    def _source_contains_docstring(src, doc):
        """
        Check whether the source *src* contains the docstring *doc*.

        This is is helper function to skip displaying the docstring if the
        source already contains it, avoiding repetition of information.
        """
        try:
            (def_node,) = ast.parse(dedent(src)).body
            return ast.get_docstring(def_node) == doc  # type: ignore[arg-type]
        except Exception:
            # The source can become invalid or even non-existent (because it
            # is re-fetched from the source file) so the above code fail in
            # arbitrary ways.
            return False

    def psearch(self,pattern,ns_table,ns_search=[],
                ignore_case=False,show_all=False, *, list_types=False):
        """Search namespaces with wildcards for objects.

        Arguments:

        - pattern: string containing shell-like wildcards to use in namespace
          searches and optionally a type specification to narrow the search to
          objects of that type.

        - ns_table: dict of name->namespaces for search.

        Optional arguments:

          - ns_search: list of namespace names to include in search.

          - ignore_case(False): make the search case-insensitive.

          - show_all(False): show all names, including those starting with
            underscores.

          - list_types(False): list all available object types for object matching.
        """
        # print('ps pattern:<%r>' % pattern)  # dbg

        # defaults
        type_pattern = 'all'
        filter = ''

        # list all object types
        if list_types:
            page.page('\n'.join(sorted(typestr2type)))
            return

        cmds = pattern.split()
        len_cmds  =  len(cmds)
        if len_cmds == 1:
            # Only filter pattern given
            filter = cmds[0]
        elif len_cmds == 2:
            # Both filter and type specified
            filter,type_pattern = cmds
        else:
            raise ValueError('invalid argument string for psearch: <%s>' %
                             pattern)

        # filter search namespaces
        for name in ns_search:
            if name not in ns_table:
                raise ValueError('invalid namespace <%s>. Valid names: %s' %
                                 (name,ns_table.keys()))

        # print('type_pattern:',type_pattern)  # dbg
        search_result, namespaces_seen = set(), set()
        for ns_name in ns_search:
            ns = ns_table[ns_name]
            # Normally, locals and globals are the same, so we just check one.
            if id(ns) in namespaces_seen:
                continue
            namespaces_seen.add(id(ns))
            tmp_res = list_namespace(ns, type_pattern, filter,
                                    ignore_case=ignore_case, show_all=show_all)
            search_result.update(tmp_res)

        page.page('\n'.join(sorted(search_result)))


def _render_signature(obj_signature, obj_name) -> str:
    """
    This was mostly taken from inspect.Signature.__str__.
    Look there for the comments.
    The only change is to add linebreaks when this gets too long.
    """
    result = []
    pos_only = False
    kw_only = True
    for param in obj_signature.parameters.values():
        if param.kind == inspect.Parameter.POSITIONAL_ONLY:
            pos_only = True
        elif pos_only:
            result.append('/')
            pos_only = False

        if param.kind == inspect.Parameter.VAR_POSITIONAL:
            kw_only = False
        elif param.kind == inspect.Parameter.KEYWORD_ONLY and kw_only:
            result.append('*')
            kw_only = False

        result.append(str(param))

    if pos_only:
        result.append('/')

    # add up name, parameters, braces (2), and commas
    if len(obj_name) + sum(len(r) + 2 for r in result) > 75:
        # This doesn’t fit behind “Signature: ” in an inspect window.
        rendered = '{}(\n{})'.format(obj_name, ''.join(
            '    {},\n'.format(r) for r in result)
        )
    else:
        rendered = '{}({})'.format(obj_name, ', '.join(result))

    if obj_signature.return_annotation is not inspect._empty:
        anno = inspect.formatannotation(obj_signature.return_annotation)
        rendered += ' -> {}'.format(anno)

    return rendered
