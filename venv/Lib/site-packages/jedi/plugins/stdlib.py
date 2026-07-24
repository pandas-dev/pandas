"""
Implementations of standard library functions, because it's not possible to
understand them with Jedi.

To add a new implementation, create a function and add it to the
``_implemented`` dict at the bottom of this module.

Note that this module exists only to implement very specific functionality in
the standard library. The usual way to understand the standard library is the
compiled module that returns the types for C-builtins.
"""
import parso
import os

from jedi import debug
from jedi.inference.utils import safe_property
from jedi.inference.helpers import get_str_or_none
from jedi.inference.arguments import iterate_argument_clinic, ParamIssue, \
    repack_with_argument_clinic, AbstractArguments, TreeArgumentsWrapper
from jedi.inference import analysis
from jedi.inference import compiled
from jedi.inference.value.instance import \
    AnonymousMethodExecutionContext, MethodExecutionContext
from jedi.inference.base_value import ContextualizedNode, \
    NO_VALUES, ValueSet, ValueWrapper, LazyValueWrapper
from jedi.inference.value import ClassValue, ModuleValue
from jedi.inference.value.decorator import Decoratee
from jedi.inference.value.klass import (
    DataclassWrapper,
    DataclassDecorator,
    DataclassTransformer,
)
from jedi.inference.value.function import FunctionMixin
from jedi.inference.value import iterable
from jedi.inference.lazy_value import LazyTreeValue, LazyKnownValue, \
    LazyKnownValues
from jedi.inference.names import ValueName
from jedi.inference.filters import AttributeOverwrite, publish_method, \
    ParserTreeFilter, DictFilter
from jedi.inference.signature import SignatureWrapper


# Copied from Python 3.6's stdlib.
_NAMEDTUPLE_CLASS_TEMPLATE = """\
_property = property
_tuple = tuple
from operator import itemgetter as _itemgetter
from collections import OrderedDict

class {typename}(tuple):
    __slots__ = ()

    _fields = {field_names!r}

    def __new__(_cls, {arg_list}):
        'Create new instance of {typename}({arg_list})'
        return _tuple.__new__(_cls, ({arg_list}))

    @classmethod
    def _make(cls, iterable, new=tuple.__new__, len=len):
        'Make a new {typename} object from a sequence or iterable'
        result = new(cls, iterable)
        if len(result) != {num_fields:d}:
            raise TypeError('Expected {num_fields:d} arguments, got %d' % len(result))
        return result

    def _replace(_self, **kwds):
        'Return a new {typename} object replacing specified fields with new values'
        result = _self._make(map(kwds.pop, {field_names!r}, _self))
        if kwds:
            raise ValueError('Got unexpected field names: %r' % list(kwds))
        return result

    def __repr__(self):
        'Return a nicely formatted representation string'
        return self.__class__.__name__ + '({repr_fmt})' % self

    def _asdict(self):
        'Return a new OrderedDict which maps field names to their values.'
        return OrderedDict(zip(self._fields, self))

    def __getnewargs__(self):
        'Return self as a plain tuple.  Used by copy and pickle.'
        return tuple(self)

    # These methods were added by Jedi.
    # __new__ doesn't really work with Jedi. So adding this to nametuples seems
    # like the easiest way.
    def __init__(self, {arg_list}):
        'A helper function for namedtuple.'
        self.__iterable = ({arg_list})

    def __iter__(self):
        for i in self.__iterable:
            yield i

    def __getitem__(self, y):
        return self.__iterable[y]

{field_defs}
"""

_NAMEDTUPLE_FIELD_TEMPLATE = '''\
    {name} = _property(_itemgetter({index:d}), doc='Alias for field number {index:d}')
'''


def execute(callback):
    def wrapper(value, arguments):
        def call():
            return callback(value, arguments=arguments)

        try:
            obj_name = value.name.string_name
        except AttributeError:
            pass
        else:
            p = value.parent_context
            if p is not None and p.is_builtins_module():
                module_name = 'builtins'
            elif p is not None and p.is_module():
                module_name = p.py__name__()
            else:
                return call()

            if value.is_bound_method() or value.is_instance():
                # value can be an instance for example if it is a partial
                # object.
                return call()

            # for now we just support builtin functions.
            try:
                func = _implemented[module_name][obj_name]
            except KeyError:
                pass
            else:
                return func(value, arguments=arguments, callback=call)  # type: ignore
        return call()

    return wrapper


def _follow_param(inference_state, arguments, index):
    try:
        key, lazy_value = list(arguments.unpack())[index]
    except IndexError:
        return NO_VALUES
    else:
        return lazy_value.infer()


def argument_clinic(clinic_string, want_value=False, want_context=False,
                    want_arguments=False, want_inference_state=False,
                    want_callback=False):
    """
    Works like Argument Clinic (PEP 436), to validate function params.
    """

    def f(func):
        def wrapper(value, arguments, callback):
            try:
                args = tuple(iterate_argument_clinic(
                    value.inference_state, arguments, clinic_string))
            except ParamIssue:
                return NO_VALUES

            debug.dbg('builtin start %s' % value, color='MAGENTA')
            kwargs = {}
            if want_context:
                kwargs['context'] = arguments.context
            if want_value:
                kwargs['value'] = value
            if want_inference_state:
                kwargs['inference_state'] = value.inference_state
            if want_arguments:
                kwargs['arguments'] = arguments
            if want_callback:
                kwargs['callback'] = callback
            result = func(*args, **kwargs)
            debug.dbg('builtin end: %s', result, color='MAGENTA')
            return result

        return wrapper
    return f


@argument_clinic('iterator[, default], /', want_inference_state=True)
def builtins_next(iterators, defaults, inference_state):
    # TODO theoretically we have to check here if something is an iterator.
    # That is probably done by checking if it's not a class.
    return defaults | iterators.py__getattribute__('__next__').execute_with_values()


@argument_clinic('iterator[, default], /')
def builtins_iter(iterators_or_callables, defaults):
    # TODO implement this if it's a callable.
    return iterators_or_callables.py__getattribute__('__iter__').execute_with_values()


@argument_clinic('object, name[, default], /')
def builtins_getattr(objects, names, defaults=None):
    # follow the first param
    for value in objects:
        for name in names:
            string = get_str_or_none(name)
            if string is None:
                debug.warning('getattr called without str')
                continue
            else:
                return value.py__getattribute__(string)
    return NO_VALUES


@argument_clinic('object[, bases, dict], /')
def builtins_type(objects, bases, dicts):
    if bases or dicts:
        # It's a type creation... maybe someday...
        return NO_VALUES
    else:
        return objects.py__class__()


class SuperInstance(LazyValueWrapper):
    """To be used like the object ``super`` returns."""
    def __init__(self, inference_state, instance):
        self.inference_state = inference_state
        self._instance = instance  # Corresponds to super().__self__

    def _get_bases(self):
        return self._instance.py__class__().py__bases__()

    def _get_wrapped_value(self):
        objs = self._get_bases()[0].infer().execute_with_values()
        if not objs:
            # This is just a fallback and will only be used, if it's not
            # possible to find a class
            return self._instance
        return next(iter(objs))

    def get_filters(self, origin_scope=None):
        for b in self._get_bases():
            for value in b.infer().execute_with_values():
                for f in value.get_filters():
                    yield f


@argument_clinic('[type[, value]], /', want_context=True)
def builtins_super(types, objects, context):
    instance = None
    if isinstance(context, AnonymousMethodExecutionContext):
        instance = context.instance
    elif isinstance(context, MethodExecutionContext):
        instance = context.instance
    if instance is None:
        return NO_VALUES
    return ValueSet({SuperInstance(instance.inference_state, instance)})


class ReversedObject(AttributeOverwrite):
    def __init__(self, reversed_obj, iter_list):
        super().__init__(reversed_obj)
        self._iter_list = iter_list

    def py__iter__(self, contextualized_node=None):
        return self._iter_list

    @publish_method('__next__')
    def _next(self, arguments):
        return ValueSet.from_sets(
            lazy_value.infer() for lazy_value in self._iter_list
        )


@argument_clinic('sequence, /', want_value=True, want_arguments=True)
def builtins_reversed(sequences, value, arguments):
    # While we could do without this variable (just by using sequences), we
    # want static analysis to work well. Therefore we need to generated the
    # values again.
    key, lazy_value = next(arguments.unpack())
    cn = None
    if isinstance(lazy_value, LazyTreeValue):
        cn = ContextualizedNode(lazy_value.context, lazy_value.data)
    ordered = list(sequences.iterate(cn))

    # Repack iterator values and then run it the normal way. This is
    # necessary, because `reversed` is a function and autocompletion
    # would fail in certain cases like `reversed(x).__iter__` if we
    # just returned the result directly.
    seq, = value.inference_state.typing_module.py__getattribute__('Iterator').execute_with_values()
    return ValueSet([ReversedObject(seq, list(reversed(ordered)))])


@argument_clinic('value, type, /', want_arguments=True, want_inference_state=True)
def builtins_isinstance(objects, types, arguments, inference_state):
    bool_results = set()
    for o in objects:
        cls = o.py__class__()
        try:
            cls.py__bases__
        except AttributeError:
            # This is temporary. Everything should have a class attribute in
            # Python?! Maybe we'll leave it here, because some numpy objects or
            # whatever might not.
            bool_results = set([True, False])
            break

        mro = list(cls.py__mro__())

        for cls_or_tup in types:
            if cls_or_tup.is_class():
                bool_results.add(cls_or_tup in mro)
            elif cls_or_tup.name.string_name == 'tuple' \
                    and cls_or_tup.get_root_context().is_builtins_module():
                # Check for tuples.
                classes = ValueSet.from_sets(
                    lazy_value.infer()
                    for lazy_value in cls_or_tup.iterate()
                )
                bool_results.add(any(cls in mro for cls in classes))
            else:
                _, lazy_value = list(arguments.unpack())[1]
                if isinstance(lazy_value, LazyTreeValue):
                    node = lazy_value.data
                    message = 'TypeError: isinstance() arg 2 must be a ' \
                              'class, type, or tuple of classes and types, ' \
                              'not %s.' % cls_or_tup
                    analysis.add(lazy_value.context, 'type-error-isinstance', node, message)

    return ValueSet(
        compiled.builtin_from_name(inference_state, str(b))
        for b in bool_results
    )


class StaticMethodObject(ValueWrapper):
    def py__get__(self, instance, class_value):
        return ValueSet([self._wrapped_value])


@argument_clinic('sequence, /')
def builtins_staticmethod(functions):
    return ValueSet(StaticMethodObject(f) for f in functions)


class ClassMethodObject(ValueWrapper):
    def __init__(self, class_method_obj, function):
        super().__init__(class_method_obj)
        self._function = function

    def py__get__(self, instance, class_value):
        return ValueSet([
            ClassMethodGet(__get__, class_value, self._function)
            for __get__ in self._wrapped_value.py__getattribute__('__get__')
        ])


class ClassMethodGet(ValueWrapper):
    def __init__(self, get_method, klass, function):
        super().__init__(get_method)
        self._class = klass
        self._function = function

    def get_signatures(self):
        return [sig.bind(self._function) for sig in self._function.get_signatures()]

    def py__call__(self, arguments):
        return self._function.execute(ClassMethodArguments(self._class, arguments))


class ClassMethodArguments(TreeArgumentsWrapper):
    def __init__(self, klass, arguments):
        super().__init__(arguments)
        self._class = klass

    def unpack(self, func=None):
        yield None, LazyKnownValue(self._class)
        for values in self._wrapped_arguments.unpack(func):
            yield values


@argument_clinic('sequence, /', want_value=True, want_arguments=True)
def builtins_classmethod(functions, value, arguments):
    return ValueSet(
        ClassMethodObject(class_method_object, function)
        for class_method_object in value.py__call__(arguments=arguments)
        for function in functions
    )


class PropertyObject(AttributeOverwrite, ValueWrapper):
    api_type = 'property'

    def __init__(self, property_obj, function):
        super().__init__(property_obj)
        self._function = function

    def py__get__(self, instance, class_value):
        if instance is None:
            return ValueSet([self])
        return self._function.execute_with_values(instance)

    @publish_method('deleter')
    @publish_method('getter')
    @publish_method('setter')
    def _return_self(self, arguments):
        return ValueSet({self})


@argument_clinic('func, /', want_callback=True)
def builtins_property(functions, callback):
    return ValueSet(
        PropertyObject(property_value, function)
        for property_value in callback()
        for function in functions
    )


def collections_namedtuple(value, arguments, callback):
    """
    Implementation of the namedtuple function.

    This has to be done by processing the namedtuple class template and
    inferring the result.

    """
    inference_state = value.inference_state

    # Process arguments
    name = 'jedi_unknown_namedtuple'
    for c in _follow_param(inference_state, arguments, 0):
        x = get_str_or_none(c)
        if x is not None:
            name = x
            break

    # TODO here we only use one of the types, we should use all.
    param_values = _follow_param(inference_state, arguments, 1)
    if not param_values:
        return NO_VALUES
    _fields = list(param_values)[0]
    string = get_str_or_none(_fields)
    if string is not None:
        fields = string.replace(',', ' ').split()
    elif isinstance(_fields, iterable.Sequence):
        fields = [
            get_str_or_none(v)
            for lazy_value in _fields.py__iter__()
            for v in lazy_value.infer()
        ]
        fields = [f for f in fields if f is not None]
    else:
        return NO_VALUES

    # Build source code
    code = _NAMEDTUPLE_CLASS_TEMPLATE.format(
        typename=name,
        field_names=tuple(fields),
        num_fields=len(fields),
        arg_list=repr(tuple(fields)).replace("'", "")[1:-1],
        repr_fmt='',
        field_defs='\n'.join(_NAMEDTUPLE_FIELD_TEMPLATE.format(index=index, name=name)
                             for index, name in enumerate(fields))
    )

    # Parse source code
    module = inference_state.grammar.parse(code)
    generated_class = next(module.iter_classdefs())
    parent_context = ModuleValue(
        inference_state, module,
        code_lines=parso.split_lines(code, keepends=True),
    ).as_context()

    return ValueSet([ClassValue(inference_state, parent_context, generated_class)])


class PartialObject(ValueWrapper):
    def __init__(self, actual_value, arguments, instance=None):
        super().__init__(actual_value)
        self._arguments = arguments
        self._instance = instance

    def _get_functions(self, unpacked_arguments):
        key, lazy_value = next(unpacked_arguments, (None, None))
        if key is not None or lazy_value is None:
            debug.warning("Partial should have a proper function %s", self._arguments)
            return None
        return lazy_value.infer()

    def get_signatures(self):
        unpacked_arguments = self._arguments.unpack()
        funcs = self._get_functions(unpacked_arguments)
        if funcs is None:
            return []

        arg_count = 0
        if self._instance is not None:
            arg_count = 1
        keys = set()
        for key, _ in unpacked_arguments:
            if key is None:
                arg_count += 1
            else:
                keys.add(key)
        return [PartialSignature(s, arg_count, keys) for s in funcs.get_signatures()]

    def py__call__(self, arguments):
        funcs = self._get_functions(self._arguments.unpack())
        if funcs is None:
            return NO_VALUES

        return funcs.execute(
            MergedPartialArguments(self._arguments, arguments, self._instance)
        )

    def py__doc__(self):
        """
        In CPython partial does not replace the docstring. However we are still
        imitating it here, because we want this docstring to be worth something
        for the user.
        """
        callables = self._get_functions(self._arguments.unpack())
        if callables is None:
            return ''
        for callable_ in callables:
            return callable_.py__doc__()
        return ''

    def py__get__(self, instance, class_value):
        return ValueSet([self])


class PartialMethodObject(PartialObject):
    def py__get__(self, instance, class_value):
        if instance is None:
            return ValueSet([self])
        return ValueSet([PartialObject(self._wrapped_value, self._arguments, instance)])


class PartialSignature(SignatureWrapper):
    def __init__(self, wrapped_signature, skipped_arg_count, skipped_arg_set):
        super().__init__(wrapped_signature)
        self._skipped_arg_count = skipped_arg_count
        self._skipped_arg_set = skipped_arg_set

    def get_param_names(self, resolve_stars=False):
        names = self._wrapped_signature.get_param_names()[self._skipped_arg_count:]
        return [n for n in names if n.string_name not in self._skipped_arg_set]


class MergedPartialArguments(AbstractArguments):
    def __init__(self, partial_arguments, call_arguments, instance=None):
        self._partial_arguments = partial_arguments
        self._call_arguments = call_arguments
        self._instance = instance

    def unpack(self, funcdef=None):
        unpacked = self._partial_arguments.unpack(funcdef)
        # Ignore this one, it's the function. It was checked before that it's
        # there.
        next(unpacked, None)
        if self._instance is not None:
            yield None, LazyKnownValue(self._instance)
        for key_lazy_value in unpacked:
            yield key_lazy_value
        for key_lazy_value in self._call_arguments.unpack(funcdef):
            yield key_lazy_value


def functools_partial(value, arguments, callback):
    return ValueSet(
        PartialObject(instance, arguments)
        for instance in value.py__call__(arguments)
    )


def functools_partialmethod(value, arguments, callback):
    return ValueSet(
        PartialMethodObject(instance, arguments)
        for instance in value.py__call__(arguments)
    )


@argument_clinic('first, /')
def _return_first_param(firsts):
    return firsts


@argument_clinic('seq')
def _random_choice(sequences):
    return ValueSet.from_sets(
        lazy_value.infer()
        for sequence in sequences
        for lazy_value in sequence.py__iter__()
    )


def _dataclass(value, arguments, callback):
    """
    Decorator entry points for dataclass.

    1. dataclass decorator declaration with parameters
    2. dataclass semantics on a class from a dataclass(-like) decorator
    """
    for c in _follow_param(value.inference_state, arguments, 0):
        if c.is_class():
            # Declare dataclass semantics on a class from a dataclass decorator
            should_generate_init = (
                # Customized decorator, init may be disabled
                value.init_param_mode
                if isinstance(value, DataclassDecorator)
                # Bare dataclass decorator, always with init mode
                else True
            )
            return ValueSet([DataclassWrapper(c, should_generate_init)])
        else:
            # @dataclass(init=False)
            # dataclass decorator customization
            return ValueSet(
                [
                    DataclassDecorator(
                        value,
                        arguments=arguments,
                        default_init=True,
                    )
                ]
            )

    return NO_VALUES


def _dataclass_transform(value, arguments, callback):
    """
    Decorator entry points for dataclass_transform.

    1. dataclass-like decorator instantiation from a dataclass_transform decorator
    2. dataclass_transform decorator declaration with parameters
    3. dataclass-like decorator declaration with parameters
    4. dataclass-like semantics on a class from a dataclass-like decorator
    """
    for c in _follow_param(value.inference_state, arguments, 0):
        if c.is_class():
            is_dataclass_transform = (
                value.name.string_name == "dataclass_transform"
                # The decorator function from dataclass_transform acting as the
                # dataclass decorator.
                and not isinstance(value, Decoratee)
                # The decorator function from dataclass_transform acting as the
                # dataclass decorator with customized parameters
                and not isinstance(value, DataclassDecorator)
            )

            if is_dataclass_transform:
                # Declare base class
                return ValueSet([DataclassTransformer(c)])
            else:
                # Declare dataclass-like semantics on a class from a
                # dataclass-like decorator
                should_generate_init = value.init_param_mode
                return ValueSet([DataclassWrapper(c, should_generate_init)])
        elif c.is_function():
            # dataclass-like decorator instantiation:
            # @dataclass_transform
            # def create_model()
            return ValueSet(
                [
                    DataclassDecorator(
                        value,
                        arguments=arguments,
                        default_init=True,
                    )
                ]
            )
        elif (
            # @dataclass_transform
            # def create_model(): pass
            # @create_model(init=...)
            isinstance(value, Decoratee)
        ):
            # dataclass (or like) decorator customization
            return ValueSet(
                [
                    DataclassDecorator(
                        value,
                        arguments=arguments,
                        default_init=value._wrapped_value.init_param_mode,
                    )
                ]
            )
        else:
            # dataclass_transform decorator with parameters; nothing impactful
            return ValueSet([value])
    return NO_VALUES


class ItemGetterCallable(ValueWrapper):
    def __init__(self, instance, args_value_set):
        super().__init__(instance)
        self._args_value_set = args_value_set

    @repack_with_argument_clinic('item, /')
    def py__call__(self, item_value_set):
        value_set = NO_VALUES
        for args_value in self._args_value_set:
            lazy_values = list(args_value.py__iter__())
            if len(lazy_values) == 1:
                # TODO we need to add the contextualized value.
                value_set |= item_value_set.get_item(lazy_values[0].infer(), None)
            else:
                value_set |= ValueSet([iterable.FakeList(
                    self._wrapped_value.inference_state,
                    [
                        LazyKnownValues(item_value_set.get_item(lazy_value.infer(), None))
                        for lazy_value in lazy_values
                    ],
                )])
        return value_set


@argument_clinic('func, /')
def _functools_wraps(funcs):
    return ValueSet(WrapsCallable(func) for func in funcs)


class WrapsCallable(ValueWrapper):
    # XXX this is not the correct wrapped value, it should be a weird
    #     partials object, but it doesn't matter, because it's always used as a
    #     decorator anyway.
    @repack_with_argument_clinic('func, /')
    def py__call__(self, funcs):
        return ValueSet({Wrapped(func, self._wrapped_value) for func in funcs})


class Wrapped(ValueWrapper, FunctionMixin):
    def __init__(self, func, original_function):
        super().__init__(func)
        self._original_function = original_function

    @property
    def name(self):
        return self._original_function.name

    def get_signature_functions(self):
        return [self]


@argument_clinic('*args, /', want_value=True, want_arguments=True)
def _operator_itemgetter(args_value_set, value, arguments):
    return ValueSet([
        ItemGetterCallable(instance, args_value_set)
        for instance in value.py__call__(arguments)
    ])


def _create_string_input_function(func):
    @argument_clinic('string, /', want_value=True, want_arguments=True)
    def wrapper(strings, value, arguments):
        def iterate():
            for value in strings:
                s = get_str_or_none(value)
                if s is not None:
                    s = func(s)
                    yield compiled.create_simple_object(value.inference_state, s)
        values = ValueSet(iterate())
        if values:
            return values
        return value.py__call__(arguments)
    return wrapper


@argument_clinic('*args, /', want_callback=True)
def _os_path_join(args_set, callback):
    if len(args_set) == 1:
        string = ''
        sequence, = args_set
        is_first = True
        for lazy_value in sequence.py__iter__():
            string_values = lazy_value.infer()
            if len(string_values) != 1:
                break
            s = get_str_or_none(next(iter(string_values)))
            if s is None:
                break
            if not is_first:
                string += os.path.sep
            string += s
            is_first = False
        else:
            return ValueSet([compiled.create_simple_object(sequence.inference_state, string)])
    return callback()


_path_overrides = {
    'dirname': _create_string_input_function(os.path.dirname),
    'abspath': _create_string_input_function(os.path.abspath),
    'relpath': _create_string_input_function(os.path.relpath),
    'join': _os_path_join,
}

_implemented = {
    'builtins': {
        'getattr': builtins_getattr,
        'type': builtins_type,
        'super': builtins_super,
        'reversed': builtins_reversed,
        'isinstance': builtins_isinstance,
        'next': builtins_next,
        'iter': builtins_iter,
        'staticmethod': builtins_staticmethod,
        'classmethod': builtins_classmethod,
        'property': builtins_property,
    },
    'copy': {
        'copy': _return_first_param,
        'deepcopy': _return_first_param,
    },
    'json': {
        'load': lambda value, arguments, callback: NO_VALUES,
        'loads': lambda value, arguments, callback: NO_VALUES,
    },
    'collections': {
        'namedtuple': collections_namedtuple,
    },
    'functools': {
        'partial': functools_partial,
        'partialmethod': functools_partialmethod,
        'wraps': _functools_wraps,
    },
    '_weakref': {
        'proxy': _return_first_param,
    },
    'random': {
        'choice': _random_choice,
    },
    'operator': {
        'itemgetter': _operator_itemgetter,
    },
    'abc': {
        # Not sure if this is necessary, but it's used a lot in typeshed and
        # it's for now easier to just pass the function.
        'abstractmethod': _return_first_param,
    },
    'typing': {
        # The _alias function just leads to some annoying type inference.
        # Therefore, just make it return nothing, which leads to the stubs
        # being used instead. This only matters for 3.7+.
        '_alias': lambda value, arguments, callback: NO_VALUES,
        # runtime_checkable doesn't really change anything and is just
        # adding logs for infering stuff, so we can safely ignore it.
        'runtime_checkable': lambda value, arguments, callback: NO_VALUES,
        # Python 3.11+
        'dataclass_transform': _dataclass_transform,
    },
    'typing_extensions': {
        # Python <3.11
        'dataclass_transform': _dataclass_transform,
    },
    'dataclasses': {
        # For now this works at least better than Jedi trying to understand it.
        'dataclass': _dataclass
    },
    'posixpath': _path_overrides,
    'ntpath': _path_overrides,
}


def get_metaclass_filters(func):
    def wrapper(cls, metaclasses, is_instance):
        for metaclass in metaclasses:
            if metaclass.py__name__() == 'EnumMeta' \
                    and metaclass.get_root_context().py__name__() == 'enum':
                filter_ = ParserTreeFilter(parent_context=cls.as_context())
                return [DictFilter({
                    name.string_name: EnumInstance(cls, name).name
                    for name in filter_.values()
                })]
        return func(cls, metaclasses, is_instance)
    return wrapper


class EnumInstance(LazyValueWrapper):
    def __init__(self, cls, name):
        self.inference_state = cls.inference_state
        self._cls = cls  # Corresponds to super().__self__
        self._name = name
        self.tree_node = self._name.tree_name

    @safe_property
    def name(self):
        return ValueName(self, self._name.tree_name)

    def _get_wrapped_value(self):
        n = self._name.string_name
        if n.startswith('__') and n.endswith('__') or self._name.api_type == 'function':
            inferred = self._name.infer()
            if inferred:
                return next(iter(inferred))
            o, = self.inference_state.builtins_module.py__getattribute__('object')
            return o

        value, = self._cls.execute_with_values()
        return value

    def get_filters(self, origin_scope=None):
        yield DictFilter(dict(
            name=compiled.create_simple_object(self.inference_state, self._name.string_name).name,
            value=self._name,
        ))
        for f in self._get_wrapped_value().get_filters():
            yield f


# Make sure tuple[...] behaves like Tuple[...]
class TupleClassWrapper(ValueWrapper):
    def py__getitem__(self, index_value_set, contextualized_node):
        return self.inference_state.typing_tuple().py__getitem__(
            index_value_set,
            contextualized_node,
        )


# Make sure type[...] behaves like Type[...]
class TypeClassWrapper(ValueWrapper):
    def py__getitem__(self, index_value_set, contextualized_node):
        return self.inference_state.typing_type().py__getitem__(
            index_value_set,
            contextualized_node,
        )


def tree_name_to_values(func):
    def wrapper(inference_state, context, tree_name):
        if tree_name.value == 'sep' \
                and context.is_module() and context.py__name__() in ('posixpath', 'ntpath'):
            return ValueSet({
                compiled.create_simple_object(inference_state, os.path.sep),
            })
        if tree_name.value == 'tuple' \
                and context.is_module() and context.py__name__() == 'builtins':
            tup, = func(inference_state, context, tree_name)
            return ValueSet([TupleClassWrapper(tup)])
        if tree_name.value == 'type' \
                and context.is_module() and context.py__name__() == 'builtins':
            tup, = func(inference_state, context, tree_name)
            return ValueSet([TypeClassWrapper(tup)])
        return func(inference_state, context, tree_name)
    return wrapper
