"""
Like described in the :mod:`parso.python.tree` module,
there's a need for an ast like module to represent the states of parsed
modules.

But now there are also structures in Python that need a little bit more than
that. An ``Instance`` for example is only a ``Class`` before it is
instantiated. This class represents these cases.

So, why is there also a ``Class`` class here? Well, there are decorators and
they change classes in Python 3.

Representation modules also define "magic methods". Those methods look like
``py__foo__`` and are typically mappable to the Python equivalents ``__call__``
and others. Here's a list:

====================================== ========================================
**Method**                             **Description**
-------------------------------------- ----------------------------------------
py__call__(arguments: Array)           On callable objects, returns types.
py__bool__()                           Returns True/False/None; None means that
                                       there's no certainty.
py__bases__()                          Returns a list of base classes.
py__iter__()                           Returns a generator of a set of types.
py__class__()                          Returns the class of an instance.
py__simple_getitem__(index: int/str)   Returns a a set of types of the index.
                                       Can raise an IndexError/KeyError.
py__getitem__(indexes: ValueSet)       Returns a a set of types of the index.
py__file__()                           Only on modules. Returns None if does
                                       not exist.
py__package__() -> List[str]           Only on modules. For the import system.
py__path__()                           Only on modules. For the import system.
py__get__(call_object)                 Only on instances. Simulates
                                       descriptors.
py__doc__()                            Returns the docstring for a value.
====================================== ========================================

"""
from __future__ import annotations

from typing import List, Optional, Tuple, TYPE_CHECKING, Any

from jedi import debug
from jedi.parser_utils import get_cached_parent_scope, expr_is_dotted, \
    function_is_property
from jedi.inference.cache import inference_state_method_cache, CachedMetaClass, \
    inference_state_method_generator_cache
from jedi.inference import compiled
from jedi.inference.lazy_value import LazyKnownValues, LazyTreeValue
from jedi.inference.filters import ParserTreeFilter
from jedi.inference.names import TreeNameDefinition, ValueName
from jedi.inference.arguments import unpack_arglist, ValuesArguments
from jedi.inference.base_value import ValueSet, iterator_to_value_set, \
    NO_VALUES, ValueWrapper
from jedi.inference.context import ClassContext
from jedi.inference.value.function import FunctionAndClassBase, FunctionMixin
from jedi.inference.value.decorator import Decoratee
from jedi.inference.gradual.generics import LazyGenericManager, TupleGenericManager
from jedi.plugins import plugin_manager
from inspect import Parameter
from jedi.inference.names import BaseTreeParamName
from jedi.inference.signature import AbstractSignature

if TYPE_CHECKING:
    from jedi.inference import InferenceState


class ClassName(TreeNameDefinition):
    def __init__(self, class_value, tree_name, name_context, apply_decorators):
        super().__init__(name_context, tree_name)
        self._apply_decorators = apply_decorators
        self._class_value = class_value

    @iterator_to_value_set
    def infer(self):
        # We're using a different value to infer, so we cannot call super().
        from jedi.inference.syntax_tree import tree_name_to_values
        inferred = tree_name_to_values(
            self.parent_context.inference_state, self.parent_context, self.tree_name)

        for result_value in inferred:
            if self._apply_decorators:
                yield from result_value.py__get__(instance=None, class_value=self._class_value)
            else:
                yield result_value

    @property
    def api_type(self):
        type_ = super().api_type
        if type_ == 'function':
            definition = self.tree_name.get_definition()
            if definition is None:
                return type_
            if function_is_property(definition):
                # This essentially checks if there is an @property before
                # the function. @property could be something different, but
                # any programmer that redefines property as something that
                # is not really a property anymore, should be shot. (i.e.
                # this is a heuristic).
                return 'property'
        return type_


class ClassFilter(ParserTreeFilter):
    def __init__(self, class_value, node_context=None, until_position=None,
                 origin_scope=None, is_instance=False):
        super().__init__(
            class_value.as_context(), node_context,
            until_position=until_position,
            origin_scope=origin_scope,
        )
        self._class_value = class_value
        self._is_instance = is_instance

    def _convert_names(self, names):
        return [
            ClassName(
                class_value=self._class_value,
                tree_name=name,
                name_context=self._node_context,
                apply_decorators=not self._is_instance,
            ) for name in names
        ]

    def _equals_origin_scope(self):
        node = self._origin_scope
        while node is not None:
            if node == self._parser_scope or node == self.parent_context:
                return True
            node = get_cached_parent_scope(self._parso_cache_node, node)
        return False

    def _access_possible(self, name):
        # Filter for name mangling of private variables like __foo
        return not name.value.startswith('__') or name.value.endswith('__') \
            or self._equals_origin_scope()

    def _filter(self, names):
        names = super()._filter(names)
        return [name for name in names if self._access_possible(name)]


def init_param_value(arg_nodes) -> Optional[bool]:
    """
    Returns:

    - ``True`` if ``@dataclass(init=True)``
    - ``False`` if ``@dataclass(init=False)``
    - ``None`` if not specified ``@dataclass()``
    """
    for arg_node in arg_nodes:
        if (
            arg_node.type == "argument"
            and arg_node.children[0].value == "init"
        ):
            if arg_node.children[2].value == "False":
                return False
            elif arg_node.children[2].value == "True":
                return True

    return None


def get_dataclass_param_names(cls) -> List[DataclassParamName]:
    """
    ``cls`` is a :class:`ClassMixin`. The type is only documented as mypy would
    complain that some fields are missing.

    .. code:: python

        @dataclass
        class A:
            a: int
            b: str = "toto"

    For the previous example, the param names would be ``a`` and ``b``.
    """
    param_names = []
    filter_ = cls.as_context().get_global_filter()
    for name in sorted(filter_.values(), key=lambda name: name.start_pos):
        d = name.tree_name.get_definition()
        annassign = d.children[1]
        if d.type == 'expr_stmt' and annassign.type == 'annassign':
            node = annassign.children[1]
            if node.type == "atom_expr" and node.children[0].value == "ClassVar":
                continue

            if len(annassign.children) < 4:
                default = None
            else:
                default = annassign.children[3]

            param_names.append(DataclassParamName(
                parent_context=cls.parent_context,
                tree_name=name.tree_name,
                annotation_node=annassign.children[1],
                default_node=default,
            ))
    return param_names


class ClassMixin:
    tree_node: Any
    parent_context: Any
    inference_state: InferenceState
    py__bases__: Any
    get_metaclasses: Any
    get_metaclass_filters: Any
    get_metaclass_signatures: Any
    list_type_vars: Any

    def is_class(self):
        return True

    def is_class_mixin(self):
        return True

    def py__call__(self, arguments):
        from jedi.inference.value import TreeInstance

        from jedi.inference.gradual.typing import TypedDict
        if self.is_typeddict():
            return ValueSet([TypedDict(self)])
        return ValueSet([TreeInstance(self.inference_state, self.parent_context, self, arguments)])

    def py__class__(self):
        return compiled.builtin_from_name(self.inference_state, 'type')

    @property
    def name(self):
        return ValueName(self, self.tree_node.name)

    def py__name__(self):
        return self.name.string_name

    @inference_state_method_generator_cache()
    def py__mro__(self):
        mro = [self]
        yield self
        # TODO Do a proper mro resolution. Currently we are just listing
        # classes. However, it's a complicated algorithm.
        for lazy_cls in self.py__bases__():
            # TODO there's multiple different mro paths possible if this yields
            # multiple possibilities. Could be changed to be more correct.
            for cls in lazy_cls.infer():
                # TODO detect for TypeError: duplicate base class str,
                # e.g.  `class X(str, str): pass`
                try:
                    mro_method = cls.py__mro__
                except AttributeError:
                    # TODO add a TypeError like:
                    """
                    >>> class Y(lambda: test): pass
                    Traceback (most recent call last):
                      File "<stdin>", line 1, in <module>
                    TypeError: function() argument 1 must be code, not str
                    >>> class Y(1): pass
                    Traceback (most recent call last):
                      File "<stdin>", line 1, in <module>
                    TypeError: int() takes at most 2 arguments (3 given)
                    """
                    debug.warning('Super class of %s is not a class: %s', self, cls)
                else:
                    for cls_new in mro_method():
                        if cls_new not in mro:
                            mro.append(cls_new)
                            yield cls_new

    def get_filters(self, origin_scope=None, is_instance=False,
                    include_metaclasses=True, include_type_when_class=True):
        if include_metaclasses:
            metaclasses = self.get_metaclasses()
            if metaclasses:
                yield from self.get_metaclass_filters(metaclasses, is_instance)

        for cls in self.py__mro__():
            if cls.is_compiled():
                yield from cls.get_filters(is_instance=is_instance)
            else:
                yield ClassFilter(
                    self, node_context=cls.as_context(),
                    origin_scope=origin_scope,
                    is_instance=is_instance
                )
        if not is_instance and include_type_when_class:
            from jedi.inference.compiled import builtin_from_name
            type_ = builtin_from_name(self.inference_state, 'type')
            if type_ != self:
                # We are not using execute_with_values here, because the
                # plugin function for type would get executed instead of an
                # instance creation.
                args = ValuesArguments([])
                for instance in type_.py__call__(args):
                    instance_filters = instance.get_filters()
                    # Filter out self filters
                    next(instance_filters, None)
                    next(instance_filters, None)
                    x = next(instance_filters, None)
                    assert x is not None
                    yield x

    def _has_dataclass_transform_metaclasses(self) -> Tuple[bool, Optional[bool]]:
        for meta in self.get_metaclasses():  # type: ignore[attr-defined]
            if (
                isinstance(meta, Decoratee)
                # Internal leakage :|
                and isinstance(meta._wrapped_value, DataclassTransformer)
            ):
                return True, meta._wrapped_value.init_mode_from_new()

        return False, None

    def _get_dataclass_transform_signatures(self) -> List[DataclassSignature]:
        """
        Returns: A non-empty list if the class has dataclass semantics else an
        empty list.

        The dataclass-like semantics will be assumed for any class that directly
        or indirectly derives from the decorated class or uses the decorated
        class as a metaclass.
        """
        param_names = []
        is_dataclass_transform = False
        default_init_mode: Optional[bool] = None
        for cls in reversed(list(self.py__mro__())):
            if not is_dataclass_transform:

                # If dataclass_transform is applied to a class, dataclass-like semantics
                # will be assumed for any class that directly or indirectly derives from
                # the decorated class or uses the decorated class as a metaclass.
                if (
                    isinstance(cls, DataclassTransformer)
                    and cls.init_mode_from_init_subclass
                ):
                    is_dataclass_transform = True
                    default_init_mode = cls.init_mode_from_init_subclass

                elif (
                    # Some object like CompiledValues would not be compatible
                    isinstance(cls, ClassMixin)
                ):
                    is_dataclass_transform, default_init_mode = (
                        cls._has_dataclass_transform_metaclasses()
                    )

                # Attributes on the decorated class and its base classes are not
                # considered to be fields.
                if is_dataclass_transform:
                    continue

            # All inherited classes behave like dataclass semantics
            if (
                is_dataclass_transform
                and isinstance(cls, ClassValue)
                and (
                    cls.init_param_mode()
                    or (cls.init_param_mode() is None and default_init_mode)
                )
            ):
                param_names.extend(
                    get_dataclass_param_names(cls)
                )

        if is_dataclass_transform:
            return [DataclassSignature(cls, param_names)]
        else:
            return []

    def get_signatures(self):
        # Since calling staticmethod without a function is illegal, the Jedi
        # plugin doesn't return anything. Therefore call directly and get what
        # we want: An instance of staticmethod.
        metaclasses = self.get_metaclasses()
        if metaclasses:
            sigs = self.get_metaclass_signatures(metaclasses)
            if sigs:
                return sigs
        args = ValuesArguments([])
        instance = self.py__call__(args)
        init_funcs = init_or_new_func(instance)

        dataclass_sigs = self._get_dataclass_transform_signatures()
        if dataclass_sigs:
            return dataclass_sigs
        else:
            return [sig.bind(self) for sig in init_funcs.get_signatures()]

    def _as_context(self):
        return ClassContext(self)

    def get_type_hint(self, add_class_info=True):
        if add_class_info:
            return 'Type[%s]' % self.py__name__()
        return self.py__name__()

    @inference_state_method_cache(default=False)
    def is_typeddict(self):
        # TODO Do a proper mro resolution. Currently we are just listing
        # classes. However, it's a complicated algorithm.
        from jedi.inference.gradual.typing import TypedDictClass
        for lazy_cls in self.py__bases__():
            if not isinstance(lazy_cls, LazyTreeValue):
                return False
            tree_node = lazy_cls.data
            # Only resolve simple classes, stuff like Iterable[str] are more
            # intensive to resolve and if generics are involved, we know it's
            # not a TypedDict.
            if not expr_is_dotted(tree_node):
                return False

            for cls in lazy_cls.infer():
                if isinstance(cls, TypedDictClass):
                    return True
                try:
                    method = cls.is_typeddict
                except AttributeError:
                    # We're only dealing with simple classes, so just returning
                    # here should be fine. This only happens with e.g. compiled
                    # classes.
                    return False
                else:
                    if method():
                        return True
        return False

    def py__getitem__(self, index_value_set, contextualized_node):
        from jedi.inference.gradual.base import GenericClass
        if not index_value_set:
            debug.warning('Class indexes inferred to nothing. Returning class instead')
            return ValueSet([self])
        return ValueSet(
            GenericClass(
                self,
                LazyGenericManager(
                    context_of_index=contextualized_node.context,
                    index_value=index_value,
                )
            )
            for index_value in index_value_set
        )

    def with_generics(self, generics_tuple):
        from jedi.inference.gradual.base import GenericClass
        return GenericClass(
            self,
            TupleGenericManager(generics_tuple)
        )

    def define_generics(self, type_var_dict):
        from jedi.inference.gradual.base import GenericClass

        def remap_type_vars():
            """
            The TypeVars in the resulting classes have sometimes different names
            and we need to check for that, e.g. a signature can be:

            def iter(iterable: Iterable[_T]) -> Iterator[_T]: ...

            However, the iterator is defined as Iterator[_T_co], which means it has
            a different type var name.
            """
            for type_var in self.list_type_vars():
                yield type_var_dict.get(type_var.py__name__(), NO_VALUES)

        if type_var_dict:
            return ValueSet([GenericClass(
                self,
                TupleGenericManager(tuple(remap_type_vars()))
            )])
        return ValueSet({self})


def init_or_new_func(value):
    init_funcs = value.py__getattribute__('__init__')
    if len(init_funcs) == 1:
        init = next(iter(init_funcs))
        try:
            class_context = init.class_context
        except AttributeError:
            pass
        else:
            # In the case where we are on object.__init__, we try to use
            # __new__.
            if class_context.get_root_context().is_builtins_module() \
                    and init.class_context.name.string_name == "object":
                return value.py__getattribute__('__new__')
    return init_funcs


class DataclassParamName(BaseTreeParamName):
    """
    Represent a field declaration on a class with dataclass semantics.
    """

    def __init__(self, parent_context, tree_name, annotation_node, default_node):
        super().__init__(parent_context, tree_name)
        self.annotation_node = annotation_node
        self.default_node = default_node

    def get_kind(self):
        return Parameter.POSITIONAL_OR_KEYWORD

    def infer(self):
        if self.annotation_node is None:
            return NO_VALUES
        else:
            return self.parent_context.infer_node(self.annotation_node)


class DataclassSignature(AbstractSignature):
    """
    It represents the ``__init__`` signature of a class with dataclass semantics.

    .. code:: python

    """
    def __init__(self, value, param_names):
        super().__init__(value)
        self._param_names = param_names

    def get_param_names(self, resolve_stars=False):
        return self._param_names


class DataclassDecorator(ValueWrapper, FunctionMixin):
    """
    A dataclass(-like) decorator with custom parameters.

    .. code:: python

        @dataclass(init=True) # this
        class A: ...

        @dataclass_transform
        def create_model(*, init=False): pass

        @create_model(init=False) # or this
        class B: ...
    """

    def __init__(self, function, arguments, default_init: bool = True):
        """
        Args:
            function: Decoratee | function
            arguments: The parameters to the dataclass function decorator
            default_init: Boolean to indicate the default init value
        """
        super().__init__(function)
        argument_init = self._init_param_value(arguments)
        self.init_param_mode = (
            argument_init if argument_init is not None else default_init
        )

    def _init_param_value(self, arguments) -> Optional[bool]:
        if not arguments.argument_node:
            return None

        arg_nodes = (
            arguments.argument_node.children
            if arguments.argument_node.type == "arglist"
            else [arguments.argument_node]
        )

        return init_param_value(arg_nodes)


class DataclassTransformer(ValueWrapper, ClassMixin):
    """
    A class decorated with the ``dataclass_transform`` decorator. dataclass-like
    semantics will be assumed for any class that directly or indirectly derives
    from the decorated class or uses the decorated class as a metaclass.
    Attributes on the decorated class and its base classes are not considered to
    be fields.
    """
    def __init__(self, wrapped_value):
        super().__init__(wrapped_value)

    def init_mode_from_new(self) -> bool:
        """Default value if missing is ``True``"""
        new_methods = self._wrapped_value.py__getattribute__("__new__")

        if not new_methods:
            return True

        new_method = list(new_methods)[0]

        for param in new_method.get_param_names():
            if (
                param.string_name == "init"
                and param.default_node
                and param.default_node.type == "keyword"
            ):
                if param.default_node.value == "False":
                    return False
                elif param.default_node.value == "True":
                    return True

        return True

    @property
    def init_mode_from_init_subclass(self) -> Optional[bool]:
        # def __init_subclass__(cls) -> None: ... is hardcoded in the typeshed
        # so the extra parameters can not be inferred.
        return True


class DataclassWrapper(ValueWrapper, ClassMixin):
    """
    A class with dataclass semantics from a decorator. The init parameters are
    only from the current class and parent classes decorated where the ``init``
    parameter was ``True``.

    .. code:: python

        @dataclass
        class A: ... # this

        @dataclass_transform
        def create_model(): pass

        @create_model()
        class B: ... # or this
    """

    def __init__(
        self, wrapped_value, should_generate_init: bool
    ):
        super().__init__(wrapped_value)
        self.should_generate_init = should_generate_init

    def get_signatures(self):
        param_names = []
        for cls in reversed(list(self.py__mro__())):
            if (
                isinstance(cls, DataclassWrapper)
                and cls.should_generate_init
            ):
                param_names.extend(get_dataclass_param_names(cls))
        return [DataclassSignature(cls, param_names)]


class ClassValue(ClassMixin, FunctionAndClassBase, metaclass=CachedMetaClass):
    api_type = 'class'

    @inference_state_method_cache()
    def list_type_vars(self):
        found = []
        arglist = self.tree_node.get_super_arglist()
        if arglist is None:
            return []

        for stars, node in unpack_arglist(arglist):
            if stars:
                continue  # These are not relevant for this search.

            from jedi.inference.gradual.annotation import find_unknown_type_vars
            for type_var in find_unknown_type_vars(self.parent_context, node):
                if type_var not in found:
                    # The order matters and it's therefore a list.
                    found.append(type_var)
        return found

    def _get_bases_arguments(self):
        arglist = self.tree_node.get_super_arglist()
        if arglist:
            from jedi.inference import arguments
            return arguments.TreeArguments(self.inference_state, self.parent_context, arglist)
        return None

    @inference_state_method_cache(default=())
    def py__bases__(self):
        args = self._get_bases_arguments()
        if args is not None:
            lst = [value for key, value in args.unpack() if key is None]
            if lst:
                return lst

        if self.py__name__() == 'object' \
                and self.parent_context.is_builtins_module():
            return []
        return [LazyKnownValues(
            self.inference_state.builtins_module.py__getattribute__('object')
        )]

    @plugin_manager.decorate()
    def get_metaclass_filters(self, metaclasses, is_instance):
        debug.warning('Unprocessed metaclass %s', metaclasses)
        return []

    @inference_state_method_cache(default=NO_VALUES)
    def get_metaclasses(self):
        args = self._get_bases_arguments()
        if args is not None:
            m = [value for key, value in args.unpack() if key == 'metaclass']
            metaclasses = ValueSet.from_sets(lazy_value.infer() for lazy_value in m)
            metaclasses = ValueSet(m for m in metaclasses if m.is_class())
            if metaclasses:
                return metaclasses

        for lazy_base in self.py__bases__():
            for value in lazy_base.infer():
                if value.is_class():
                    values = value.get_metaclasses()
                    if values:
                        return values
        return NO_VALUES

    def init_param_mode(self) -> Optional[bool]:
        """
        It returns ``True`` if ``class X(init=False):`` else ``False``.
        """
        bases_arguments = self._get_bases_arguments()

        if bases_arguments is None:
            return None

        if bases_arguments.argument_node.type != "arglist":
            # If it is not inheriting from the base model and having
            # extra parameters, then init behavior is not changed.
            return None

        return init_param_value(bases_arguments.argument_node.children)

    @plugin_manager.decorate()
    def get_metaclass_signatures(self, metaclasses):
        return []
