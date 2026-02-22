from __future__ import annotations

import operator
import re
from enum import Enum
from typing import TYPE_CHECKING, Literal, NewType, TypeVar

from sphinx.errors import PycodeError
from sphinx.events import EventManager
from sphinx.ext.autodoc._directive_options import _AutoDocumenterOptions
from sphinx.ext.autodoc._dynamic._loader import _load_object_by_name
from sphinx.ext.autodoc._dynamic._mock import ismock, undecorate
from sphinx.ext.autodoc._property_types import _ClassDefProperties, _ModuleProperties
from sphinx.ext.autodoc._sentinels import ALL, INSTANCE_ATTR, SLOTS_ATTR
from sphinx.ext.autodoc._shared import LOGGER
from sphinx.locale import __
from sphinx.pycode import ModuleAnalyzer
from sphinx.util import inspect
from sphinx.util.docstrings import separate_metadata
from sphinx.util.inspect import (
    getannotations,
    getdoc,
    getmro,
    getslots,
    isclass,
    isenumclass,
    safe_getattr,
    unwrap_all,
)
from sphinx.util.typing import AnyTypeAliasType

if TYPE_CHECKING:
    from collections.abc import Iterable, Iterator, Mapping, MutableSet, Sequence, Set
    from typing import Any, Literal

    from sphinx.environment import _CurrentDocument
    from sphinx.events import EventManager
    from sphinx.ext.autodoc._directive_options import _AutoDocumenterOptions
    from sphinx.ext.autodoc._property_types import _AutodocObjType, _ItemProperties
    from sphinx.ext.autodoc._sentinels import (
        ALL_T,
        EMPTY_T,
        INSTANCE_ATTR_T,
        SLOTS_ATTR_T,
    )
    from sphinx.ext.autodoc._shared import _AttrGetter, _AutodocConfig

special_member_re = re.compile(r'^__\S+__$')


class ObjectMember:
    """A member of object.

    This is used for the result of `_get_members_to_document()` to
    represent each member of the object.
    """

    __slots__ = '__name__', 'object', 'docstring', 'class_'

    __name__: str
    object: Any
    docstring: Sequence[str] | None
    class_: Any
    skipped: bool

    def __init__(
        self,
        name: str,
        obj: INSTANCE_ATTR_T | SLOTS_ATTR_T | Any,
        *,
        docstring: Sequence[str] | None = None,
        class_: Any = None,
    ) -> None:
        self.__name__ = name
        self.object = obj
        self.docstring = docstring
        self.class_ = class_

    def __repr__(self) -> str:
        return (
            f'ObjectMember('
            f'name={self.__name__!r}, '
            f'obj={self.object!r}, '
            f'docstring={self.docstring!r}, '
            f'class_={self.class_!r}'
            f')'
        )


def _gather_members(
    *,
    want_all: bool,
    indent: str,
    analyzer_order: dict[str, int],
    attr_docs: dict[tuple[str, str], list[str]],
    config: _AutodocConfig,
    current_document: _CurrentDocument,
    events: EventManager,
    get_attr: _AttrGetter,
    options: _AutoDocumenterOptions,
    parent_modname: str,
    props: _ItemProperties,
    ref_context: Mapping[str, str | None],
    reread_always: MutableSet[str],
) -> list[tuple[_ItemProperties, bool, str]]:
    """Generate reST for member documentation.

    If *want_all* is True, document all members, else those given by
    *self.options.members*.
    """
    if props.obj_type not in {'module', 'class', 'exception'}:
        msg = 'must be implemented in subclasses'
        raise NotImplementedError(msg)
    assert isinstance(props, (_ModuleProperties, _ClassDefProperties))

    indent += '   ' * (props.obj_type != 'module')

    # set current namespace for finding members
    current_document.autodoc_module = props.module_name
    if props.parts:
        current_document.autodoc_class = props.parts[0]

    inherited_members = frozenset(options.inherited_members or ())
    found_members = _get_members_to_document(
        want_all=want_all,
        get_attr=get_attr,
        class_signature=config.autodoc_class_signature,
        inherit_docstrings=config.autodoc_inherit_docstrings,
        props=props,
        opt_members=options.members or (),
        inherited_members=inherited_members,
        opt_private_members=options.private_members,
        opt_special_members=options.special_members,
        ignore_module_all=bool(options.ignore_module_all),
        attr_docs=attr_docs,
    )
    filtered_members = _filter_members(
        found_members,
        want_all=want_all,
        events=events,
        get_attr=get_attr,
        class_signature=config.autodoc_class_signature,
        inherit_docstrings=config.autodoc_inherit_docstrings,
        options=options,
        props=props,
        inherited_members=inherited_members,
        exclude_members=options.exclude_members,
        special_members=options.special_members,
        private_members=options.private_members,
        undoc_members=options.undoc_members,
        attr_docs=attr_docs,
    )
    # document non-skipped members
    member_documenters: list[tuple[_ItemProperties, bool, str]] = []
    for member_name, member, is_attr in filtered_members:
        # prefer the object type with the highest priority
        obj_type = _best_object_type_for_member(
            member=member,
            member_name=member_name,
            is_attr=is_attr,
            parent_obj_type=props.obj_type,
            parent_props=props,
        )
        if not obj_type:
            # don't know how to document this member
            continue
        # give explicitly separated module name, so that members
        # of inner classes can be documented
        dotted_parts = '.'.join((*props.parts, member_name))
        full_name = f'{props.module_name}::{dotted_parts}'

        # We now try to import all objects before ordering them. This is to
        # avoid possible circular imports if we were to import objects after
        # their associated documenters have been sorted.
        member_props = _load_object_by_name(
            name=full_name,
            objtype=obj_type,
            current_document=current_document,
            config=config,
            events=events,
            get_attr=get_attr,
            options=options,
            parent_modname=parent_modname,
            ref_context=ref_context,
            reread_always=reread_always,
        )
        if member_props is None:
            continue
        member_documenters.append((member_props, is_attr, indent))

    member_order = options.member_order or config.autodoc_member_order
    member_documenters = _sort_members(
        member_documenters,
        member_order,
        ignore_module_all=bool(options.ignore_module_all),
        analyzer_order=analyzer_order,
        props=props,
    )

    # reset current objects
    current_document.autodoc_module = ''
    current_document.autodoc_class = ''

    return member_documenters


def _get_members_to_document(
    *,
    want_all: bool,
    get_attr: _AttrGetter,
    class_signature: Literal['mixed', 'separated'],
    inherit_docstrings: bool,
    props: _ModuleProperties | _ClassDefProperties,
    opt_members: ALL_T | Sequence[str],
    inherited_members: Set[str],
    opt_private_members: ALL_T | Sequence[str] | None,
    opt_special_members: ALL_T | Sequence[str] | None,
    ignore_module_all: bool,
    attr_docs: dict[tuple[str, str], list[str]],
) -> list[ObjectMember]:
    """Find out which members are documentable

    If *want_all* is True, return all members.  Else, only return those
    members given by *self.options.members* (which may also be None).

    Filter the given member list.

    Members are skipped if

    - they are private (except if given explicitly or the private-members
      option is set)
    - they are special methods (except if given explicitly or the
      special-members option is set)
    - they are undocumented (except if the undoc-members option is set)

    The user can override the skipping decision by connecting to the
    ``autodoc-skip-member`` event.
    """
    wanted_members: ALL_T | Set[str]
    if want_all:
        if (
            props.obj_type == 'module'
            and not ignore_module_all
            and props.all is not None
        ):
            wanted_members = frozenset(props.all)
        else:
            wanted_members = ALL
    else:
        # specific members given
        assert opt_members is not ALL

        # Merge :private-members: and :special-members: into :members:
        combined_members = set(opt_members)
        if opt_private_members is not None and opt_private_members is not ALL:
            combined_members.update(opt_private_members)
        if opt_special_members is not None and opt_special_members is not ALL:
            combined_members.update(opt_special_members)
        if class_signature == 'separated' and props.obj_type in {'class', 'exception'}:
            combined_members |= {'__new__', '__init__'}  # show __init__() method
        wanted_members = frozenset(combined_members)

    object_members_map: dict[str, ObjectMember] = {}
    if props.obj_type == 'module':
        for name in dir(props._obj):
            try:
                value = safe_getattr(props._obj, name, None)
                if ismock(value):
                    value = undecorate(value)
                if name in wanted_members:
                    object_members_map[name] = ObjectMember(
                        name, value, docstring=attr_docs.get(('', name), [])
                    )
            except AttributeError:
                continue

        # annotation only member (e.g. attr: int)
        for name in inspect.getannotations(props._obj):
            if name not in object_members_map and name in wanted_members:
                object_members_map[name] = ObjectMember(
                    name, INSTANCE_ATTR, docstring=attr_docs.get(('', name), [])
                )

        obj_members_seq = list(object_members_map.values())
    elif props.obj_type in {'class', 'exception'}:
        # the members directly defined in the class
        obj_dict = get_attr(props._obj, '__dict__', {})

        # enum members
        if isenumclass(props._obj):
            for name, defining_class, value in _filter_enum_dict(
                props._obj, get_attr, obj_dict
            ):
                # the order of occurrence of *name* matches obj's MRO,
                # allowing inherited attributes to be shadowed correctly
                if unmangled := unmangle(defining_class, name):
                    if unmangled in wanted_members:
                        object_members_map[unmangled] = ObjectMember(
                            unmangled, value, class_=defining_class
                        )

        # members in __slots__
        try:
            subject___slots__ = getslots(props._obj)
            if subject___slots__:
                for name, subject_docstring in subject___slots__.items():
                    if name not in wanted_members:
                        continue
                    if isinstance(subject_docstring, str):
                        subject_doclines = subject_docstring.splitlines()
                    else:
                        subject_doclines = None
                    object_members_map[name] = ObjectMember(
                        name,
                        SLOTS_ATTR,
                        class_=props._obj,
                        docstring=subject_doclines,
                    )
        except (TypeError, ValueError):
            pass

        # other members
        for name in dir(props._obj):
            try:
                value = get_attr(props._obj, name)
                if ismock(value):
                    value = undecorate(value)

                unmangled = unmangle(props._obj, name)
                if (
                    unmangled
                    and unmangled not in object_members_map
                    and unmangled in wanted_members
                ):
                    if name in obj_dict:
                        object_members_map[unmangled] = ObjectMember(
                            unmangled, value, class_=props._obj
                        )
                    else:
                        object_members_map[unmangled] = ObjectMember(unmangled, value)
            except AttributeError:
                continue

        try:
            for cls in getmro(props._obj):
                try:
                    modname = safe_getattr(cls, '__module__')
                    qualname = safe_getattr(cls, '__qualname__')
                except AttributeError:
                    qualname = None
                    analyzer = None
                else:
                    try:
                        analyzer = ModuleAnalyzer.for_module(modname)
                        analyzer.analyze()
                    except PycodeError:
                        analyzer = None

                # annotation only member (ex. attr: int)
                for name in getannotations(cls):
                    unmangled = unmangle(cls, name)
                    if (
                        unmangled
                        and unmangled not in object_members_map
                        and unmangled in wanted_members
                    ):
                        if analyzer and (qualname, unmangled) in analyzer.attr_docs:
                            attr_docstring = analyzer.attr_docs[qualname, unmangled]
                        else:
                            attr_docstring = None
                        object_members_map[unmangled] = ObjectMember(
                            unmangled,
                            INSTANCE_ATTR,
                            class_=cls,
                            docstring=attr_docstring,
                        )

                # append or complete instance attributes (cf. self.attr1) if analyzer knows
                if analyzer:
                    for (ns, name), attr_docstring in analyzer.attr_docs.items():
                        if ns == qualname and name not in object_members_map:
                            # otherwise unknown instance attribute
                            if name in wanted_members:
                                object_members_map[name] = ObjectMember(
                                    name,
                                    INSTANCE_ATTR,
                                    class_=cls,
                                    docstring=attr_docstring,
                                )
                        elif (
                            ns == qualname
                            and attr_docstring
                            and not object_members_map[name].docstring
                        ):
                            if cls != props._obj and not inherit_docstrings:
                                # If we are in the MRO of the class and not the class itself,
                                # and we do not want to inherit docstrings, then skip setting
                                # the docstring below
                                continue
                            # attribute is already known, because dir(props._obj)
                            # enumerates it. But it has no docstring yet
                            object_members_map[name].docstring = attr_docstring
        except AttributeError:
            pass

        if want_all and not inherited_members:
            obj_members_seq = [
                m for m in object_members_map.values() if m.class_ == props._obj
            ]
        else:
            obj_members_seq = list(object_members_map.values())
    else:
        raise ValueError

    if not want_all and opt_members is not ALL:
        for name in opt_members:
            if name in object_members_map:
                continue
            msg = __(
                'attribute %s is listed in :members: but is missing '
                'as it was not found in object %r'
            )
            LOGGER.warning(msg, name, props._obj, type='autodoc')
    return obj_members_seq


def _filter_members(
    obj_members_seq: Iterable[ObjectMember],
    *,
    want_all: bool,
    events: EventManager,
    get_attr: _AttrGetter,
    options: _AutoDocumenterOptions,
    props: _ModuleProperties | _ClassDefProperties,
    class_signature: Literal['mixed', 'separated'],
    inherit_docstrings: bool,
    inherited_members: Set[str],
    exclude_members: EMPTY_T | Set[str] | None,
    special_members: ALL_T | Sequence[str] | None,
    private_members: ALL_T | Sequence[str] | None,
    undoc_members: Literal[True] | None,
    attr_docs: dict[tuple[str, str], list[str]],
) -> Iterator[tuple[str, Any, bool]]:
    # search for members in source code too
    namespace = props.dotted_parts  # will be empty for modules

    # process members and determine which to skip
    for obj in obj_members_seq:
        member_name = obj.__name__
        member_obj = obj.object
        has_attr_doc = (namespace, member_name) in attr_docs
        try:
            keep = _should_keep_member(
                member_name=member_name,
                member_obj=member_obj,
                member_docstring=obj.docstring,
                member_cls=obj.class_,
                get_attr=get_attr,
                has_attr_doc=has_attr_doc,
                class_signature=class_signature,
                inherit_docstrings=inherit_docstrings,
                inherited_members=inherited_members,
                parent=props._obj,
                want_all=want_all,
                exclude_members=exclude_members,
                special_members=special_members,
                private_members=private_members,
                undoc_members=undoc_members,
            )
        except Exception as exc:
            LOGGER.warning(
                __(
                    'autodoc: failed to determine %s.%s (%r) to be documented, '
                    'the following exception was raised:\n%s'
                ),
                props.full_name,
                member_name,
                member_obj,
                exc,
                type='autodoc',
            )
            keep = False

        # give the user a chance to decide whether this member
        # should be skipped
        if events is not None:
            # let extensions preprocess docstrings
            skip_member = events.emit_firstresult(
                'autodoc-skip-member',
                props.obj_type,
                member_name,
                member_obj,
                not keep,
                options,
            )
            if skip_member is not None:
                keep = not skip_member

        if keep:
            # if is_attr is True, the member is documented as an attribute
            is_attr = member_obj is INSTANCE_ATTR or has_attr_doc
            yield member_name, member_obj, is_attr


def _best_object_type_for_member(
    member: Any,
    member_name: str,
    is_attr: bool,
    *,
    parent_obj_type: str,
    parent_props: _ItemProperties | None,
) -> _AutodocObjType | None:
    """Return the best object type that supports documenting *member*."""
    filtered = []

    # Don't document submodules automatically: 'module' is never returned.

    try:
        if isinstance(member, type) and issubclass(member, BaseException):
            # priority must be higher than 'class'
            filtered.append((20, 'exception'))
    except TypeError as exc:
        # It's possible for a member to be considered a type, but fail
        # issubclass checks due to not being a class. For example:
        # https://github.com/sphinx-doc/sphinx/issues/11654#issuecomment-1696790436
        msg = f'Failed to discern if member {member} is a BaseException subclass.'
        raise ValueError(msg) from exc

    if isinstance(member, type) or (is_attr and isinstance(member, (NewType, TypeVar))):
        # priority must be higher than 'function', 'class', and 'attribute'
        # as NewType can be an attribute and is a class.
        filtered.append((15, 'class'))

    if parent_obj_type in {'class', 'exception'}:
        if inspect.isproperty(member):
            # priority must be higher than 'attribute'
            filtered.append((11, 'property'))

        # See _get_documenter() in autosummary, parent_props might be None.
        elif parent_props is not None:
            # Support for class properties. Note: these only work on Python 3.9.
            __dict__ = safe_getattr(parent_props._obj, '__dict__', {})
            obj = __dict__.get(member_name)
            if isinstance(obj, classmethod) and inspect.isproperty(obj.__func__):
                # priority must be higher than 'attribute'
                filtered.append((11, 'property'))

    if parent_obj_type != 'module':
        if inspect.isattributedescriptor(member) or not (
            inspect.isroutine(member) or isinstance(member, type)
        ):
            # priority must be higher than 'method', else it will recognise
            # some non-data descriptors as methods
            filtered.append((10, 'attribute'))

    if inspect.isroutine(member) and parent_obj_type != 'module':
        # priority must be higher than 'function'
        filtered.append((1, 'method'))

    if (
        inspect.isfunction(member)
        or inspect.isbuiltin(member)
        or (inspect.isroutine(member) and parent_obj_type == 'module')
    ):
        # supports functions, builtins and bound methods exported
        # at the module level
        filtered.extend(((0, 'function'), (-1, 'decorator')))

    if isinstance(member, AnyTypeAliasType):
        filtered.append((0, 'type'))

    if parent_obj_type == 'module' and is_attr:
        filtered.append((-10, 'data'))

    if filtered:
        # return the highest priority object type
        return max(filtered, key=operator.itemgetter(0))[1]  # type: ignore[return-value]
    return None


def _sort_members(
    documenters: list[tuple[_ItemProperties, bool, str]],
    order: Literal['alphabetical', 'bysource', 'groupwise'],
    *,
    ignore_module_all: bool,
    analyzer_order: dict[str, int],
    props: _ItemProperties,
) -> list[tuple[_ItemProperties, bool, str]]:
    """Sort the given member list."""
    if order == 'groupwise':
        # sort by group; alphabetically within groups
        def group_order(entry: tuple[_ItemProperties, bool, str]) -> tuple[int, str]:
            return entry[0]._groupwise_order_key, entry[0].full_name

        documenters.sort(key=group_order)
    elif order == 'bysource':
        if (
            isinstance(props, _ModuleProperties)
            and not ignore_module_all
            and (module_all := props.all)
        ):
            # Sort by __all__
            module_all_idx = {name: idx for idx, name in enumerate(module_all)}
            module_all_len = len(module_all)

            def source_order(entry: tuple[_ItemProperties, bool, str]) -> int:
                fullname = entry[0].dotted_parts
                return module_all_idx.get(fullname, module_all_len)

            documenters.sort(key=source_order)

        # By default, member discovery order matches source order,
        # as dicts are insertion-ordered from Python 3.7.
        elif analyzer_order:
            # sort by source order, by virtue of the module analyzer
            order_len = len(analyzer_order)

            def source_order(entry: tuple[_ItemProperties, bool, str]) -> int:
                fullname = entry[0].dotted_parts
                return analyzer_order.get(fullname, order_len)

            documenters.sort(key=source_order)
    else:  # alphabetical
        documenters.sort(key=lambda entry: entry[0].full_name)

    return documenters


def unmangle(subject: Any, name: str) -> str | None:
    """Unmangle the given name."""
    try:
        if isclass(subject) and not name.endswith('__'):
            prefix = f'_{subject.__name__}__'
            if name.startswith(prefix):
                return name.replace(prefix, '__', 1)
            else:
                for cls in subject.__mro__:
                    prefix = f'_{cls.__name__}__'
                    if name.startswith(prefix):
                        # mangled attribute defined in parent class
                        return None
    except AttributeError:
        pass

    return name


def _filter_enum_dict(
    enum_class: type[Enum],
    attrgetter: _AttrGetter,
    enum_class_dict: Mapping[str, object],
) -> Iterator[tuple[str, type, Any]]:
    """Find the attributes to document of an enumeration class.

    The output consists of triplets ``(attribute name, defining class, value)``
    where the attribute name can appear more than once during the iteration
    but with different defining class. The order of occurrence is guided by
    the MRO of *enum_class*.
    """
    # attributes that were found on a mixin type or the data type
    candidate_in_mro: set[str] = set()
    # sunder names that were picked up (and thereby allowed to be redefined)
    # see: https://docs.python.org/3/howto/enum.html#supported-dunder-names
    sunder_names = {
        '_name_',
        '_value_',
        '_missing_',
        '_order_',
        '_generate_next_value_',
    }
    # attributes that can be picked up on a mixin type or the enum's data type
    public_names = {'name', 'value', *object.__dict__, *sunder_names}
    # names that are ignored by default
    ignore_names = Enum.__dict__.keys() - public_names

    def should_ignore(name: str, value: Any) -> bool:
        if name in sunder_names:
            return _is_native_enum_api(value, name)
        return name in ignore_names

    sentinel = object()

    def query(name: str, defining_class: type) -> tuple[str, type, Any] | None:
        value = attrgetter(enum_class, name, sentinel)
        if value is not sentinel:
            return name, defining_class, value
        return None

    # attributes defined on a parent type, possibly shadowed later by
    # the attributes defined directly inside the enumeration class
    for parent in enum_class.__mro__:
        if parent in {enum_class, Enum, object}:
            continue

        parent_dict = attrgetter(parent, '__dict__', {})
        for name, value in parent_dict.items():
            if should_ignore(name, value):
                continue

            candidate_in_mro.add(name)
            if (item := query(name, parent)) is not None:
                yield item

    # exclude members coming from the native Enum unless
    # they were redefined on a mixin type or the data type
    excluded_members = Enum.__dict__.keys() - candidate_in_mro
    yield from filter(
        None,
        (
            query(name, enum_class)
            for name in enum_class_dict
            if name not in excluded_members
        ),
    )

    # check if allowed members from ``Enum`` were redefined at the enum level
    special_names = sunder_names | public_names
    special_names &= enum_class_dict.keys()
    special_names &= Enum.__dict__.keys()
    for name in special_names:
        if (
            not _is_native_enum_api(enum_class_dict[name], name)
            and (item := query(name, enum_class)) is not None
        ):
            yield item


def _is_native_enum_api(obj: object, name: str) -> bool:
    """Check whether *obj* is the same as ``Enum.__dict__[name]``."""
    return unwrap_all(obj) is unwrap_all(Enum.__dict__[name])


def _should_keep_member(
    *,
    member_name: str,
    member_obj: Any,
    member_docstring: Sequence[str] | None,
    member_cls: Any,
    get_attr: _AttrGetter,
    has_attr_doc: bool,
    class_signature: Literal['mixed', 'separated'],
    inherit_docstrings: bool,
    inherited_members: Set[str],
    parent: Any,
    want_all: bool,
    exclude_members: EMPTY_T | Set[str] | None,
    special_members: ALL_T | Sequence[str] | None,
    private_members: ALL_T | Sequence[str] | None,
    undoc_members: Literal[True] | None,
) -> bool:
    if member_docstring:
        # hack for ClassDocumenter to inject docstring
        doclines: Sequence[str] | None = member_docstring
    else:
        doc = getdoc(
            member_obj,
            get_attr,
            inherit_docstrings,
            parent,
            member_name,
        )
        # Ignore non-string __doc__
        doclines = doc.splitlines() if isinstance(doc, str) else None

        # if the member __doc__ is the same as self's __doc__, it's just
        # inherited and therefore not the member's doc
        cls = get_attr(member_obj, '__class__', None)
        if cls:
            cls_doc = get_attr(cls, '__doc__', None)
            if cls_doc == doc:
                doclines = None

    if doclines is not None:
        doc, metadata = separate_metadata('\n'.join(doclines))
    else:
        doc = ''
        metadata = {}
    has_doc = bool(doc or undoc_members)

    if 'private' in metadata:
        # consider a member private if docstring has "private" metadata
        is_private = True
    elif 'public' in metadata:
        # consider a member public if docstring has "public" metadata
        is_private = False
    else:
        is_private = member_name.startswith('_')

    if ismock(member_obj) and not has_attr_doc:
        # mocked module or object
        return False

    if exclude_members and member_name in exclude_members:
        # remove members given by exclude-members
        return False

    if not want_all:
        # keep documented attributes
        return has_doc or has_attr_doc

    is_filtered_inherited_member = _is_filtered_inherited_member(
        member_name,
        member_cls=member_cls,
        parent=parent,
        inherited_members=inherited_members,
        get_attr=get_attr,
    )

    if special_member_re.match(member_name):
        # special __methods__
        if member_name == '__doc__' or is_filtered_inherited_member:
            return False
        if special_members and member_name in special_members:
            return has_doc
        if (
            class_signature == 'separated'
            and member_name in {'__new__', '__init__'}
            and inspect.isclass(parent)
        ):
            return has_doc  # show __init__() method
        return False

    if is_private:
        if has_attr_doc or has_doc:
            if private_members is None:  # NoQA: SIM114
                return False
            elif has_doc and is_filtered_inherited_member:
                return False
            return member_name in private_members
        return False

    if has_attr_doc:
        # keep documented attributes
        return True

    if is_filtered_inherited_member:
        return False

    # ignore undocumented members if :undoc-members: is not given
    return has_doc


def _is_filtered_inherited_member(
    member_name: str,
    *,
    member_cls: Any,
    parent: Any,
    inherited_members: Set[str],
    get_attr: _AttrGetter,
) -> bool:
    if not inspect.isclass(parent):
        return False

    seen = set()
    for cls in parent.__mro__:
        if member_name in cls.__dict__:
            seen.add(cls)
        if (
            cls.__name__ in inherited_members
            and cls != parent
            and any(issubclass(potential_child, cls) for potential_child in seen)
        ):
            # given member is a member of specified *super class*
            return True
        if member_cls is cls:
            return False
        if member_name in cls.__dict__:
            return False
        if member_name in get_attr(cls, '__annotations__', {}):
            return False
    return False
