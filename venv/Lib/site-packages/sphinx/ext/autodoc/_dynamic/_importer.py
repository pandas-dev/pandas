"""Importer utilities for autodoc"""

from __future__ import annotations

import contextlib
import importlib
import os
import sys
import traceback
import typing
from importlib.abc import FileLoader
from importlib.machinery import EXTENSION_SUFFIXES
from importlib.util import decode_source, find_spec, module_from_spec, spec_from_loader
from pathlib import Path
from typing import TYPE_CHECKING

from sphinx.errors import PycodeError
from sphinx.ext.autodoc._dynamic._mock import ismock, mock, undecorate
from sphinx.ext.autodoc._sentinels import RUNTIME_INSTANCE_ATTRIBUTE, UNINITIALIZED_ATTR
from sphinx.ext.autodoc._shared import LOGGER
from sphinx.pycode import ModuleAnalyzer
from sphinx.util import inspect
from sphinx.util.inspect import isclass, safe_getattr
from sphinx.util.typing import get_type_hints

if TYPE_CHECKING:
    from collections.abc import Mapping, Sequence
    from importlib.machinery import ModuleSpec
    from types import ModuleType
    from typing import Any, Protocol

    from sphinx.ext.autodoc._property_types import _AutodocObjType

    class _AttrGetter(Protocol):
        def __call__(self, obj: Any, name: str, default: Any = ..., /) -> Any: ...


_NATIVE_SUFFIXES: frozenset[str] = frozenset({'.pyx', *EXTENSION_SUFFIXES})


class _ImportedObject:
    #: module containing the object to document
    module: ModuleType | None

    #: parent/owner of the object to document
    parent: Any

    #: name of the object to document
    object_name: str

    #: object to document
    obj: Any

    def __init__(
        self,
        *,
        module: ModuleType | None = None,
        parent: Any,
        object_name: str = '',
        obj: Any,
    ) -> None:
        self.module = module
        self.parent = parent
        self.object_name = object_name
        self.obj = obj

    def __repr__(self) -> str:
        return f'<{self.__class__.__name__} {self.__dict__}>'


def _import_object(
    *,
    get_attr: _AttrGetter = safe_getattr,
    mock_imports: Sequence[str],
    module_name: str,
    obj_path: Sequence[str],
    obj_type: _AutodocObjType,
    type_aliases: Mapping[str, str] | None,
) -> _ImportedObject | None:
    """Import the module and get the object to document."""
    try:
        with mock(mock_imports):
            im = _import_from_module_and_path(
                module_name=module_name, obj_path=obj_path, get_attr=get_attr
            )
    except ImportError as exc:
        if obj_type == 'data':
            im_ = _import_data_declaration(
                module_name=module_name,
                obj_path=obj_path,
                mock_imports=mock_imports,
                type_aliases=type_aliases,
            )
            if im_ is not None:
                return im_
        elif obj_type == 'attribute':
            im_ = _import_attribute_declaration(
                module_name=module_name,
                obj_path=obj_path,
                mock_imports=mock_imports,
                type_aliases=type_aliases,
                get_attr=get_attr,
            )
            if im_ is not None:
                return im_
        LOGGER.warning(exc.args[0], type='autodoc', subtype='import_object')
        return None

    if ismock(im.obj):
        im.obj = undecorate(im.obj)
    return im


def _import_from_module_and_path(
    *,
    module_name: str,
    obj_path: Sequence[str],
    get_attr: _AttrGetter = safe_getattr,
) -> _ImportedObject:
    obj_path = list(obj_path)
    if obj_path:
        LOGGER.debug('[autodoc] from %s import %s', module_name, '.'.join(obj_path))
    else:
        LOGGER.debug('[autodoc] import %s', module_name)

    module = None
    exc_on_importing = None
    try:
        while module is None:
            try:
                module = _import_module(module_name, try_reload=True)
                LOGGER.debug('[autodoc] import %s => %r', module_name, module)
            except ImportError as exc:
                LOGGER.debug('[autodoc] import %s => failed', module_name)
                exc_on_importing = exc
                if '.' not in module_name:
                    raise

                # retry with parent module
                module_name, _, name = module_name.rpartition('.')
                obj_path.insert(0, name)

        obj = module
        parent = None
        object_name = ''
        for attr_name in obj_path:
            parent = obj
            LOGGER.debug('[autodoc] getattr(_, %r)', attr_name)
            mangled_name = _mangle_name(obj, attr_name)
            obj = get_attr(obj, mangled_name)

            try:
                LOGGER.debug('[autodoc] => %r', obj)
            except TypeError:
                # fallback of failure on logging for broken object
                # See: https://github.com/sphinx-doc/sphinx/issues/9095
                LOGGER.debug('[autodoc] => %r', (obj,))

            object_name = attr_name
        return _ImportedObject(
            module=module,
            parent=parent,
            object_name=object_name,
            obj=obj,
        )
    except (AttributeError, ImportError) as exc:
        if isinstance(exc, AttributeError) and exc_on_importing:
            # restore ImportError
            exc = exc_on_importing

        if obj_path:
            dotted_objpath = '.'.join(obj_path)
            err_parts = [
                f'autodoc: failed to import {dotted_objpath!r} '
                f'from module {module_name!r}'
            ]
        else:
            err_parts = [f'autodoc: failed to import {module_name!r}']

        if isinstance(exc, ImportError):
            # _import_module() raises ImportError having real exception obj and
            # traceback
            real_exc = exc.args[0]
            traceback_msg = ''.join(traceback.format_exception(exc))
            if isinstance(real_exc, SystemExit):
                err_parts.append(
                    'the module executes module level statement '
                    'and it might call sys.exit().'
                )
            elif isinstance(real_exc, ImportError) and real_exc.args:
                err_parts.append(
                    f'the following exception was raised:\n{real_exc.args[0]}'
                )
            else:
                err_parts.append(
                    f'the following exception was raised:\n{traceback_msg}'
                )
        else:
            err_parts.append(
                f'the following exception was raised:\n{traceback.format_exc()}'
            )

        errmsg = '; '.join(err_parts)
        LOGGER.debug(errmsg)
        raise ImportError(errmsg) from exc


def _import_module(modname: str, try_reload: bool = False) -> Any:
    if modname in sys.modules:
        return sys.modules[modname]

    skip_pyi = bool(os.getenv('SPHINX_AUTODOC_IGNORE_NATIVE_MODULE_TYPE_STUBS', ''))
    original_module_names = frozenset(sys.modules)
    try:
        spec = find_spec(modname)
        if spec is None:
            msg = f'No module named {modname!r}'
            raise ModuleNotFoundError(msg, name=modname)  # NoQA: TRY301
        spec, pyi_path = _find_type_stub_spec(spec, modname)
        if skip_pyi or pyi_path is None:
            module = importlib.import_module(modname)
        else:
            if spec.loader is None:
                msg = 'missing loader'
                raise ImportError(msg, name=spec.name)  # NoQA: TRY301
            sys.modules[modname] = module = module_from_spec(spec)
            spec.loader.exec_module(module)
    except ImportError:
        raise
    except BaseException as exc:
        # Importing modules may cause any side effects, including
        # SystemExit, so we need to catch all errors.
        raise ImportError(exc, traceback.format_exc()) from exc
    if try_reload and os.environ.get('SPHINX_AUTODOC_RELOAD_MODULES'):
        new_modules = [m for m in sys.modules if m not in original_module_names]
        # Try reloading modules with ``typing.TYPE_CHECKING == True``.
        try:
            typing.TYPE_CHECKING = True  # type: ignore[misc]
            # Ignore failures; we've already successfully loaded these modules
            with contextlib.suppress(ImportError, KeyError):
                for m in new_modules:
                    mod_path = getattr(sys.modules[m], '__file__', '')
                    if mod_path and mod_path.endswith('.pyi'):
                        continue
                    _reload_module(sys.modules[m])
        finally:
            typing.TYPE_CHECKING = False  # type: ignore[misc]
        module = sys.modules[modname]
    return module


def _find_type_stub_spec(
    spec: ModuleSpec, modname: str
) -> tuple[ModuleSpec, Path | None]:
    """Try finding a spec for a PEP 561 '.pyi' stub file for native modules."""
    if spec.origin is None:
        return spec, None

    for suffix in _NATIVE_SUFFIXES:
        if not spec.origin.endswith(suffix):
            continue
        pyi_path = Path(spec.origin.removesuffix(suffix) + '.pyi')
        if not pyi_path.is_file():
            continue
        pyi_loader = _StubFileLoader(modname, path=str(pyi_path))
        pyi_spec = spec_from_loader(modname, loader=pyi_loader)
        if pyi_spec is not None:
            return pyi_spec, pyi_path
    return spec, None


class _StubFileLoader(FileLoader):
    """Load modules from ``.pyi`` stub files."""

    def get_source(self, fullname: str) -> str:
        path = self.get_filename(fullname)
        for suffix in _NATIVE_SUFFIXES:
            if not path.endswith(suffix):
                continue
            path = path.removesuffix(suffix) + '.pyi'
        try:
            source_bytes = self.get_data(path)
        except OSError as exc:
            raise ImportError from exc
        return decode_source(source_bytes)


def _reload_module(module: ModuleType) -> Any:
    """Call importlib.reload(module), convert exceptions to ImportError"""
    try:
        return importlib.reload(module)
    except BaseException as exc:
        # Importing modules may cause any side effects, including
        # SystemExit, so we need to catch all errors.
        raise ImportError(exc, traceback.format_exc()) from exc


def _mangle_name(subject: Any, name: str) -> str:
    """Mangle the given name."""
    try:
        if isclass(subject) and name.startswith('__') and not name.endswith('__'):
            return f'_{subject.__name__}{name}'
    except AttributeError:
        pass

    return name


def _import_data_declaration(
    *,
    module_name: str,
    obj_path: Sequence[str],
    mock_imports: Sequence[str],
    type_aliases: Mapping[str, str] | None,
) -> _ImportedObject | None:
    # annotation only instance variable (PEP-526)
    try:
        with mock(mock_imports):
            parent = _import_module(module_name)
        annotations = get_type_hints(parent, None, type_aliases, include_extras=True)
        if obj_path[-1] in annotations:
            im = _ImportedObject(
                parent=parent,
                obj=UNINITIALIZED_ATTR,
            )
            return im
    except ImportError:
        pass
    return None


def _import_attribute_declaration(
    *,
    module_name: str,
    obj_path: Sequence[str],
    mock_imports: Sequence[str],
    type_aliases: Mapping[str, str] | None,
    get_attr: _AttrGetter = safe_getattr,
) -> _ImportedObject | None:
    # Support runtime & uninitialized instance attributes.
    #
    # The former are defined in __init__() methods with doc-comments.
    # The latter are PEP-526 style annotation only annotations.
    #
    # class Foo:
    #     attr: int  #: uninitialized attribute
    #
    #     def __init__(self):
    #         self.attr = None  #: runtime attribute
    try:
        with mock(mock_imports):
            ret = _import_from_module_and_path(
                module_name=module_name, obj_path=obj_path[:-1], get_attr=get_attr
            )
        parent = ret.obj
        if _is_runtime_instance_attribute(parent=parent, obj_path=obj_path):
            im = _ImportedObject(
                parent=parent,
                obj=RUNTIME_INSTANCE_ATTRIBUTE,
            )
            return im
        elif _is_uninitialized_instance_attribute(
            parent=parent, obj_path=obj_path, type_aliases=type_aliases
        ):
            im = _ImportedObject(
                parent=parent,
                obj=UNINITIALIZED_ATTR,
            )
            return im
    except ImportError:
        pass
    return None


def _is_runtime_instance_attribute(*, parent: Any, obj_path: Sequence[str]) -> bool:
    """Check the subject is an attribute defined in __init__()."""
    # An instance variable defined in __init__().
    if _get_attribute_comment(parent=parent, obj_path=obj_path, attrname=obj_path[-1]):
        return True
    return _is_runtime_instance_attribute_not_commented(
        parent=parent, obj_path=obj_path
    )


def _is_runtime_instance_attribute_not_commented(
    *, parent: Any, obj_path: Sequence[str]
) -> bool:
    """Check the subject is an attribute defined in __init__() without comment."""
    for cls in inspect.getmro(parent):
        try:
            module = safe_getattr(cls, '__module__')
            qualname = safe_getattr(cls, '__qualname__')

            analyzer = ModuleAnalyzer.for_module(module)
            analyzer.analyze()
            if qualname and obj_path:
                key = f'{qualname}.{obj_path[-1]}'
                if key in analyzer.tagorder:
                    return True
        except (AttributeError, PycodeError):
            pass

    return False


def _get_attribute_comment(
    parent: Any, obj_path: Sequence[str], attrname: str
) -> list[str] | None:
    for cls in inspect.getmro(parent):
        try:
            module = safe_getattr(cls, '__module__')
            qualname = safe_getattr(cls, '__qualname__')

            analyzer = ModuleAnalyzer.for_module(module)
            analyzer.analyze()
            if qualname and obj_path:
                key = (qualname, attrname)
                if key in analyzer.attr_docs:
                    return list(analyzer.attr_docs[key])
        except (AttributeError, PycodeError):
            pass

    return None


def _is_uninitialized_instance_attribute(
    *, parent: Any, obj_path: Sequence[str], type_aliases: Mapping[str, str] | None
) -> bool:
    """Check the subject is an annotation only attribute."""
    annotations = get_type_hints(parent, None, type_aliases, include_extras=True)
    return obj_path[-1] in annotations
