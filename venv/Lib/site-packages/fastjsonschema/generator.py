from collections import OrderedDict
from decimal import Decimal
import re

from .exceptions import JsonSchemaValueException, JsonSchemaDefinitionException
from .indent import indent
from .ref_resolver import RefResolver


def enforce_list(variable):
    if isinstance(variable, list):
        return variable
    return [variable]


# pylint: disable=too-many-instance-attributes,too-many-public-methods
class CodeGenerator:
    """
    This class is not supposed to be used directly. Anything
    inside of this class can be changed without noticing.

    This class generates code of validation function from JSON
    schema object as string. Example:

    .. code-block:: python

        CodeGenerator(json_schema_definition).func_code
    """

    INDENT = 4  # spaces

    def __init__(self, definition, resolver=None, detailed_exceptions=True):
        self._code = []
        self._compile_regexps = {}
        self._custom_formats = {}
        self._detailed_exceptions = detailed_exceptions

        # Any extra library should be here to be imported only once.
        # Lines are imports to be printed in the file and objects
        # key-value pair to pass to compile function directly.
        self._extra_imports_lines = [
            "from decimal import Decimal",
        ]
        self._extra_imports_objects = {
            "Decimal": Decimal,
        }

        self._variables = set()
        self._indent = 0
        self._indent_last_line = None
        self._variable = None
        self._variable_name = None
        self._root_definition = definition
        self._definition = None

        # map schema URIs to validation function names for functions
        # that are not yet generated, but need to be generated
        self._needed_validation_functions = {}
        # validation function names that are already done
        self._validation_functions_done = set()

        if resolver is None:
            resolver = RefResolver.from_schema(definition, store={})
        self._resolver = resolver

        # add main function to `self._needed_validation_functions`
        self._needed_validation_functions[self._resolver.get_uri()] = self._resolver.get_scope_name()

        self._json_keywords_to_function = OrderedDict()

    @property
    def func_code(self):
        """
        Returns generated code of whole validation function as string.
        """
        self._generate_func_code()

        return '\n'.join(self._code)

    @property
    def global_state(self):
        """
        Returns global variables for generating function from ``func_code``. Includes
        compiled regular expressions and imports, so it does not have to do it every
        time when validation function is called.
        """
        self._generate_func_code()

        return dict(
            **self._extra_imports_objects,
            REGEX_PATTERNS=self._compile_regexps,
            re=re,
            JsonSchemaValueException=JsonSchemaValueException,
        )

    @property
    def global_state_code(self):
        """
        Returns global variables for generating function from ``func_code`` as code.
        Includes compiled regular expressions and imports.
        """
        self._generate_func_code()

        if not self._compile_regexps:
            return '\n'.join(self._extra_imports_lines + [
                'from fastjsonschema import JsonSchemaValueException',
                '',
                '',
            ])
        return '\n'.join(self._extra_imports_lines + [
            'import re',
            'from fastjsonschema import JsonSchemaValueException',
            '',
            '',
            'REGEX_PATTERNS = ' + serialize_regexes(self._compile_regexps),
            '',
        ])


    def _generate_func_code(self):
        if not self._code:
            self.generate_func_code()

    def generate_func_code(self):
        """
        Creates base code of validation function and calls helper
        for creating code by definition.
        """
        self.l('NoneType = type(None)')
        # Generate parts that are referenced and not yet generated
        while self._needed_validation_functions:
            # During generation of validation function, could be needed to generate
            # new one that is added again to `_needed_validation_functions`.
            # Therefore usage of while instead of for loop.
            uri, name = self._needed_validation_functions.popitem()
            self.generate_validation_function(uri, name)

    def generate_validation_function(self, uri, name):
        """
        Generate validation function for given uri with given name
        """
        self._validation_functions_done.add(uri)
        self.l('')
        with self._resolver.resolving(uri) as definition:
            with self.l('def {}(data, custom_formats={{}}, name_prefix=None):', name):
                self.generate_func_code_block(definition, 'data', 'data', clear_variables=True)
                self.l('return data')

    def generate_func_code_block(self, definition, variable, variable_name, clear_variables=False):
        """
        Creates validation rules for current definition.

        Returns the number of validation rules generated as code.
        """
        backup = self._definition, self._variable, self._variable_name
        self._definition, self._variable, self._variable_name = definition, variable, variable_name
        if clear_variables:
            backup_variables = self._variables
            self._variables = set()

        count = self._generate_func_code_block(definition)

        self._definition, self._variable, self._variable_name = backup
        if clear_variables:
            self._variables = backup_variables

        return count

    def _generate_func_code_block(self, definition):
        if not isinstance(definition, dict):
            raise JsonSchemaDefinitionException("definition must be an object")
        if '$ref' in definition:
            # needed because ref overrides any sibling keywords
            return self.generate_ref()
        else:
            return self.run_generate_functions(definition)

    def run_generate_functions(self, definition):
        """Returns the number of generate functions that were executed."""
        count = 0
        for key, func in self._json_keywords_to_function.items():
            if key in definition:
                func()
                count += 1
        return count

    def generate_ref(self):
        """
        Ref can be link to remote or local definition.

        .. code-block:: python

            {'$ref': 'http://json-schema.org/draft-04/schema#'}
            {
                'properties': {
                    'foo': {'type': 'integer'},
                    'bar': {'$ref': '#/properties/foo'}
                }
            }
        """
        with self._resolver.in_scope(self._definition['$ref']):
            name = self._resolver.get_scope_name()
            uri = self._resolver.get_uri()
            if uri not in self._validation_functions_done:
                self._needed_validation_functions[uri] = name
            # call validation function
            assert self._variable_name.startswith("data")
            path = self._variable_name[4:]
            name_arg = '(name_prefix or "data") + "{}"'.format(path)
            if '{' in name_arg:
                name_arg = name_arg + '.format(**locals())'
            self.l('{}({variable}, custom_formats, {name_arg})', name, name_arg=name_arg)


    # pylint: disable=invalid-name
    @indent
    def l(self, line, *args, **kwds):
        """
        Short-cut of line. Used for inserting line. It's formated with parameters
        ``variable``, ``variable_name`` (as ``name`` for short-cut), all keys from
        current JSON schema ``definition`` and also passed arguments in ``args``
        and named ``kwds``.

        .. code-block:: python

            self.l('if {variable} not in {enum}: raise JsonSchemaValueException("Wrong!")')

        When you want to indent block, use it as context manager. For example:

        .. code-block:: python

            with self.l('if {variable} not in {enum}:'):
                self.l('raise JsonSchemaValueException("Wrong!")')
        """
        spaces = ' ' * self.INDENT * self._indent

        name = self._variable_name
        if name:
            # Add name_prefix to the name when it is being outputted.
            assert name.startswith('data')
            name = '" + (name_prefix or "data") + "' + name[4:]
            if '{' in name:
                name = name + '".format(**locals()) + "'

        context = dict(
            self._definition if self._definition and self._definition is not True else {},
            variable=self._variable,
            name=name,
            **kwds
        )
        line = line.format(*args, **context)
        line = line.replace('\n', '\\n').replace('\r', '\\r')
        self._code.append(spaces + line)
        return line

    def e(self, string):
        """
        Short-cut of escape. Used for inserting user values into a string message.

        .. code-block:: python

            self.l('raise JsonSchemaValueException("Variable: {}")', self.e(variable))
        """
        return str(string).replace('"', '\\"')

    def exc(self, msg, *args, append_to_msg=None, rule=None):
        """
        Short-cut for creating raising exception in the code.
        """
        if not self._detailed_exceptions:
            self.l('raise JsonSchemaValueException("'+msg+'")', *args)
            return

        arg = '"'+msg+'"'
        if append_to_msg:
            arg += ' + (' + append_to_msg + ')'
        msg = 'raise JsonSchemaValueException('+arg+', value={variable}, name="{name}", definition={definition}, rule={rule})'
        definition = self._expand_refs(self._definition)
        definition_rule = self.e(definition.get(rule) if isinstance(definition, dict) else None)
        self.l(msg, *args, definition=repr(definition), rule=repr(rule), definition_rule=definition_rule)

    def _expand_refs(self, definition):
        if isinstance(definition, list):
            return [self._expand_refs(v) for v in definition]
        if not isinstance(definition, dict):
            return definition
        if "$ref" in definition and isinstance(definition["$ref"], str):
            with self._resolver.resolving(definition["$ref"]) as schema:
                return schema
        return {k: self._expand_refs(v) for k, v in definition.items()}

    def create_variable_with_length(self):
        """
        Append code for creating variable with length of that variable
        (for example length of list or dictionary) with name ``{variable}_len``.
        It can be called several times and always it's done only when that variable
        still does not exists.
        """
        variable_name = '{}_len'.format(self._variable)
        if variable_name in self._variables:
            return
        self._variables.add(variable_name)
        self.l('{variable}_len = len({variable})')

    def create_variable_keys(self):
        """
        Append code for creating variable with keys of that variable (dictionary)
        with a name ``{variable}_keys``. Similar to `create_variable_with_length`.
        """
        variable_name = '{}_keys'.format(self._variable)
        if variable_name in self._variables:
            return
        self._variables.add(variable_name)
        self.l('{variable}_keys = set({variable}.keys())')

    def create_variable_is_list(self):
        """
        Append code for creating variable with bool if it's instance of list
        with a name ``{variable}_is_list``. Similar to `create_variable_with_length`.
        """
        variable_name = '{}_is_list'.format(self._variable)
        if variable_name in self._variables:
            return
        self._variables.add(variable_name)
        self.l('{variable}_is_list = isinstance({variable}, (list, tuple))')

    def create_variable_is_dict(self):
        """
        Append code for creating variable with bool if it's instance of list
        with a name ``{variable}_is_dict``. Similar to `create_variable_with_length`.
        """
        variable_name = '{}_is_dict'.format(self._variable)
        if variable_name in self._variables:
            return
        self._variables.add(variable_name)
        self.l('{variable}_is_dict = isinstance({variable}, dict)')


def serialize_regexes(patterns_dict):
    # Unfortunately using `pprint.pformat` is causing errors
    # specially with big regexes
    regex_patterns = (
        repr(k) + ": " + repr_regex(v)
        for k, v in patterns_dict.items()
    )
    return '{\n    ' + ",\n    ".join(regex_patterns) + "\n}"


def repr_regex(regex):
    all_flags = ("A", "I", "DEBUG", "L", "M", "S", "X")
    flags = " | ".join(f"re.{f}" for f in all_flags if regex.flags & getattr(re, f))
    flags = ", " + flags if flags else ""
    return "re.compile({!r}{})".format(regex.pattern, flags)
