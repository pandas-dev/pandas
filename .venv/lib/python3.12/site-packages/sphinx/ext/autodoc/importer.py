"""Importer utilities for autodoc"""

from __future__ import annotations

import traceback
from importlib.machinery import EXTENSION_SUFFIXES
from typing import TYPE_CHECKING, NamedTuple

from sphinx.errors import PycodeError
from sphinx.ext.autodoc._dynamic._importer import (
    _find_type_stub_spec,  # NoQA: F401
    _reload_module,  # NoQA: F401
    _StubFileLoader,  # NoQA: F401
)
from sphinx.ext.autodoc._dynamic._importer import (
    _import_module as import_module,
)
from sphinx.ext.autodoc._dynamic._importer import _mangle_name as mangle
from sphinx.ext.autodoc._dynamic._member_finder import (
    _filter_enum_dict,
    _is_native_enum_api,  # NoQA: F401
    unmangle,
)
from sphinx.ext.autodoc._dynamic._mock import ismock, undecorate
from sphinx.ext.autodoc._shared import LOGGER
from sphinx.pycode import ModuleAnalyzer
from sphinx.util.inspect import (
    getannotations,
    getmro,
    getslots,
    isenumclass,
    safe_getattr,
)

if TYPE_CHECKING:
    from typing import Any, Protocol

    from sphinx.ext.autodoc._legacy_class_based._documenters import ObjectMember

    class _AttrGetter(Protocol):
        def __call__(self, obj: Any, name: str, default: Any = ..., /) -> Any: ...


_NATIVE_SUFFIXES: frozenset[str] = frozenset({'.pyx', *EXTENSION_SUFFIXES})


# Retained: legacy class-based
def import_object(
    modname: str,
    objpath: list[str],
    objtype: str = '',
    attrgetter: _AttrGetter = safe_getattr,
) -> Any:
    if objpath:
        LOGGER.debug('[autodoc] from %s import %s', modname, '.'.join(objpath))
    else:
        LOGGER.debug('[autodoc] import %s', modname)

    try:
        module = None
        exc_on_importing = None
        objpath = objpath.copy()
        while module is None:
            try:
                module = import_module(modname, try_reload=True)
                LOGGER.debug('[autodoc] import %s => %r', modname, module)
            except ImportError as exc:
                LOGGER.debug('[autodoc] import %s => failed', modname)
                exc_on_importing = exc
                if '.' in modname:
                    # retry with parent module
                    modname, name = modname.rsplit('.', 1)
                    objpath.insert(0, name)
                else:
                    raise

        obj = module
        parent = None
        object_name = None
        for attrname in objpath:
            parent = obj
            LOGGER.debug('[autodoc] getattr(_, %r)', attrname)
            mangled_name = mangle(obj, attrname)
            obj = attrgetter(obj, mangled_name)

            try:
                LOGGER.debug('[autodoc] => %r', obj)
            except TypeError:
                # fallback of failure on logging for broken object
                # See: https://github.com/sphinx-doc/sphinx/issues/9095
                LOGGER.debug('[autodoc] => %r', (obj,))

            object_name = attrname
        return [module, parent, object_name, obj]
    except (AttributeError, ImportError) as exc:
        if isinstance(exc, AttributeError) and exc_on_importing:
            # restore ImportError
            exc = exc_on_importing

        if objpath:
            errmsg = 'autodoc: failed to import %s %r from module %r' % (
                objtype,
                '.'.join(objpath),
                modname,
            )
        else:
            errmsg = f'autodoc: failed to import {objtype} {modname!r}'

        if isinstance(exc, ImportError):
            # import_module() raises ImportError having real exception obj and
            # traceback
            real_exc = exc.args[0]
            traceback_msg = traceback.format_exception(exc)
            if isinstance(real_exc, SystemExit):
                errmsg += (
                    '; the module executes module level statement '
                    'and it might call sys.exit().'
                )
            elif isinstance(real_exc, ImportError) and real_exc.args:
                errmsg += '; the following exception was raised:\n%s' % real_exc.args[0]
            else:
                errmsg += '; the following exception was raised:\n%s' % traceback_msg
        else:
            errmsg += (
                '; the following exception was raised:\n%s' % traceback.format_exc()
            )

        LOGGER.debug(errmsg)
        raise ImportError(errmsg) from exc


# Retained: legacy class-based
class Attribute(NamedTuple):
    name: str
    directly_defined: bool
    value: Any


# Retained: legacy class-based
def get_object_members(
    subject: Any,
    objpath: list[str],
    attrgetter: _AttrGetter,
    analyzer: ModuleAnalyzer | None = None,
) -> dict[str, Attribute]:
    """Get members and attributes of target object."""
    from sphinx.ext.autodoc._legacy_class_based._sentinels import INSTANCEATTR

    # the members directly defined in the class
    obj_dict = attrgetter(subject, '__dict__', {})

    members: dict[str, Attribute] = {}

    # enum members
    if isenumclass(subject):
        for name, defining_class, value in _filter_enum_dict(
            subject, attrgetter, obj_dict
        ):
            # the order of occurrence of *name* matches the subject's MRO,
            # allowing inherited attributes to be shadowed correctly
            if unmangled := unmangle(defining_class, name):
                members[unmangled] = Attribute(
                    name=unmangled,
                    directly_defined=defining_class is subject,
                    value=value,
                )

    # members in __slots__
    try:
        subject___slots__ = getslots(subject)
        if subject___slots__:
            from sphinx.ext.autodoc._legacy_class_based._sentinels import SLOTSATTR

            for name in subject___slots__:
                members[name] = Attribute(
                    name=name, directly_defined=True, value=SLOTSATTR
                )
    except (TypeError, ValueError):
        pass

    # other members
    for name in dir(subject):
        try:
            value = attrgetter(subject, name)
            directly_defined = name in obj_dict
            unmangled = unmangle(subject, name)
            if unmangled and unmangled not in members:
                members[unmangled] = Attribute(
                    name=unmangled, directly_defined=directly_defined, value=value
                )
        except AttributeError:
            continue

    # annotation only member (ex. attr: int)
    for cls in getmro(subject):
        for name in getannotations(cls):
            unmangled = unmangle(cls, name)
            if unmangled and unmangled not in members:
                members[unmangled] = Attribute(
                    name=unmangled, directly_defined=cls is subject, value=INSTANCEATTR
                )

    if analyzer:
        # append instance attributes (cf. self.attr1) if analyzer knows
        namespace = '.'.join(objpath)
        for ns, name in analyzer.find_attr_docs():
            if namespace == ns and name not in members:
                members[name] = Attribute(
                    name=name, directly_defined=True, value=INSTANCEATTR
                )

    return members


# Retained: legacy class-based
def get_class_members(
    subject: Any, objpath: Any, attrgetter: _AttrGetter, inherit_docstrings: bool = True
) -> dict[str, ObjectMember]:
    """Get members and attributes of target class."""
    from sphinx.ext.autodoc._legacy_class_based._documenters import ObjectMember
    from sphinx.ext.autodoc._legacy_class_based._sentinels import INSTANCEATTR

    # the members directly defined in the class
    obj_dict = attrgetter(subject, '__dict__', {})

    members: dict[str, ObjectMember] = {}

    # enum members
    if isenumclass(subject):
        for name, defining_class, value in _filter_enum_dict(
            subject, attrgetter, obj_dict
        ):
            # the order of occurrence of *name* matches the subject's MRO,
            # allowing inherited attributes to be shadowed correctly
            if unmangled := unmangle(defining_class, name):
                members[unmangled] = ObjectMember(
                    unmangled, value, class_=defining_class
                )

    # members in __slots__
    try:
        subject___slots__ = getslots(subject)
        if subject___slots__:
            from sphinx.ext.autodoc._legacy_class_based._sentinels import SLOTSATTR

            for name, docstring in subject___slots__.items():
                members[name] = ObjectMember(
                    name, SLOTSATTR, class_=subject, docstring=docstring
                )
    except (TypeError, ValueError):
        pass

    # other members
    for name in dir(subject):
        try:
            value = attrgetter(subject, name)
            if ismock(value):
                value = undecorate(value)

            unmangled = unmangle(subject, name)
            if unmangled and unmangled not in members:
                if name in obj_dict:
                    members[unmangled] = ObjectMember(unmangled, value, class_=subject)
                else:
                    members[unmangled] = ObjectMember(unmangled, value)
        except AttributeError:
            continue

    try:
        for cls in getmro(subject):
            try:
                modname = safe_getattr(cls, '__module__')
                qualname = safe_getattr(cls, '__qualname__')
                analyzer = ModuleAnalyzer.for_module(modname)
                analyzer.analyze()
            except AttributeError:
                qualname = None
                analyzer = None
            except PycodeError:
                analyzer = None

            # annotation only member (ex. attr: int)
            for name in getannotations(cls):
                unmangled = unmangle(cls, name)
                if unmangled and unmangled not in members:
                    if analyzer and (qualname, unmangled) in analyzer.attr_docs:
                        docstring = '\n'.join(analyzer.attr_docs[qualname, unmangled])
                    else:
                        docstring = None

                    members[unmangled] = ObjectMember(
                        unmangled, INSTANCEATTR, class_=cls, docstring=docstring
                    )

            # append or complete instance attributes (cf. self.attr1) if analyzer knows
            if analyzer:
                for (ns, name), docstring in analyzer.attr_docs.items():
                    if ns == qualname and name not in members:
                        # otherwise unknown instance attribute
                        members[name] = ObjectMember(
                            name,
                            INSTANCEATTR,
                            class_=cls,
                            docstring='\n'.join(docstring),
                        )
                    elif (
                        ns == qualname
                        and docstring
                        and isinstance(members[name], ObjectMember)
                        and not members[name].docstring
                    ):
                        if cls != subject and not inherit_docstrings:
                            # If we are in the MRO of the class and not the class itself,
                            # and we do not want to inherit docstrings, then skip setting
                            # the docstring below
                            continue
                        # attribute is already known, because dir(subject) enumerates it.
                        # But it has no docstring yet
                        members[name].docstring = '\n'.join(docstring)
    except AttributeError:
        pass

    return members
