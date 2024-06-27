import decimal
import re

from .exceptions import JsonSchemaDefinitionException
from .generator import CodeGenerator, enforce_list


JSON_TYPE_TO_PYTHON_TYPE = {
    'null': 'NoneType',
    'boolean': 'bool',
    'number': 'int, float, Decimal',
    'integer': 'int',
    'string': 'str',
    'array': 'list, tuple',
    'object': 'dict',
}

DOLLAR_FINDER = re.compile(r"(?<!\\)\$")  # Finds any un-escaped $ (including inside []-sets)


# pylint: disable=too-many-instance-attributes,too-many-public-methods
class CodeGeneratorDraft04(CodeGenerator):
    # pylint: disable=line-too-long
    # I was thinking about using ipaddress module instead of regexps for example, but it's big
    # difference in performance. With a module I got this difference: over 100 ms with a module
    # vs. 9 ms with a regex! Other modules are also ineffective or not available in standard
    # library. Some regexps are not 100% precise but good enough, fast and without dependencies.
    FORMAT_REGEXS = {
        'date-time': r'^\d{4}-[01]\d-[0-3]\d(t|T)[0-2]\d:[0-5]\d:[0-5]\d(?:\.\d+)?(?:[+-][0-2]\d:[0-5]\d|[+-][0-2]\d[0-5]\d|z|Z)\Z',
        'email': r'^[^@]+@[^@]+\.[^@]+\Z',
        'hostname': r'^(([a-zA-Z0-9]|[a-zA-Z0-9][a-zA-Z0-9\-]{0,61}[a-zA-Z0-9])\.)*([A-Za-z0-9]|[A-Za-z0-9][A-Za-z0-9\-]{0,61}[A-Za-z0-9])\Z',
        'ipv4': r'^((25[0-5]|2[0-4][0-9]|[01]?[0-9][0-9]?)\.){3}(25[0-5]|2[0-4][0-9]|[01]?[0-9][0-9]?)\Z',
        'ipv6': r'^(?:(?:[0-9A-Fa-f]{1,4}:){6}(?:[0-9A-Fa-f]{1,4}:[0-9A-Fa-f]{1,4}|(?:(?:[0-9]|[1-9][0-9]|1[0-9]{2}|2[0-4][0-9]|25[0-5])\\.){3}(?:[0-9]|[1-9][0-9]|1[0-9]{2}|2[0-4][0-9]|25[0-5]))|::(?:[0-9A-Fa-f]{1,4}:){5}(?:[0-9A-Fa-f]{1,4}:[0-9A-Fa-f]{1,4}|(?:(?:[0-9]|[1-9][0-9]|1[0-9]{2}|2[0-4][0-9]|25[0-5])\\.){3}(?:[0-9]|[1-9][0-9]|1[0-9]{2}|2[0-4][0-9]|25[0-5]))|(?:[0-9A-Fa-f]{1,4})?::(?:[0-9A-Fa-f]{1,4}:){4}(?:[0-9A-Fa-f]{1,4}:[0-9A-Fa-f]{1,4}|(?:(?:[0-9]|[1-9][0-9]|1[0-9]{2}|2[0-4][0-9]|25[0-5])\\.){3}(?:[0-9]|[1-9][0-9]|1[0-9]{2}|2[0-4][0-9]|25[0-5]))|(?:[0-9A-Fa-f]{1,4}:[0-9A-Fa-f]{1,4})?::(?:[0-9A-Fa-f]{1,4}:){3}(?:[0-9A-Fa-f]{1,4}:[0-9A-Fa-f]{1,4}|(?:(?:[0-9]|[1-9][0-9]|1[0-9]{2}|2[0-4][0-9]|25[0-5])\\.){3}(?:[0-9]|[1-9][0-9]|1[0-9]{2}|2[0-4][0-9]|25[0-5]))|(?:(?:[0-9A-Fa-f]{1,4}:){,2}[0-9A-Fa-f]{1,4})?::(?:[0-9A-Fa-f]{1,4}:){2}(?:[0-9A-Fa-f]{1,4}:[0-9A-Fa-f]{1,4}|(?:(?:[0-9]|[1-9][0-9]|1[0-9]{2}|2[0-4][0-9]|25[0-5])\\.){3}(?:[0-9]|[1-9][0-9]|1[0-9]{2}|2[0-4][0-9]|25[0-5]))|(?:(?:[0-9A-Fa-f]{1,4}:){,3}[0-9A-Fa-f]{1,4})?::[0-9A-Fa-f]{1,4}:(?:[0-9A-Fa-f]{1,4}:[0-9A-Fa-f]{1,4}|(?:(?:[0-9]|[1-9][0-9]|1[0-9]{2}|2[0-4][0-9]|25[0-5])\\.){3}(?:[0-9]|[1-9][0-9]|1[0-9]{2}|2[0-4][0-9]|25[0-5]))|(?:(?:[0-9A-Fa-f]{1,4}:){,4}[0-9A-Fa-f]{1,4})?::(?:[0-9A-Fa-f]{1,4}:[0-9A-Fa-f]{1,4}|(?:(?:[0-9]|[1-9][0-9]|1[0-9]{2}|2[0-4][0-9]|25[0-5])\\.){3}(?:[0-9]|[1-9][0-9]|1[0-9]{2}|2[0-4][0-9]|25[0-5]))|(?:(?:[0-9A-Fa-f]{1,4}:){,5}[0-9A-Fa-f]{1,4})?::[0-9A-Fa-f]{1,4}|(?:(?:[0-9A-Fa-f]{1,4}:){,6}[0-9A-Fa-f]{1,4})?::)\Z',
        'uri': r'^\w+:(\/?\/?)[^\s]+\Z',
    }

    def __init__(self, definition, resolver=None, formats={}, use_default=True, use_formats=True):
        super().__init__(definition, resolver)
        self._custom_formats = formats
        self._use_formats = use_formats
        self._use_default = use_default
        self._json_keywords_to_function.update((
            ('type', self.generate_type),
            ('enum', self.generate_enum),
            ('allOf', self.generate_all_of),
            ('anyOf', self.generate_any_of),
            ('oneOf', self.generate_one_of),
            ('not', self.generate_not),
            ('minLength', self.generate_min_length),
            ('maxLength', self.generate_max_length),
            ('pattern', self.generate_pattern),
            ('format', self.generate_format),
            ('minimum', self.generate_minimum),
            ('maximum', self.generate_maximum),
            ('multipleOf', self.generate_multiple_of),
            ('minItems', self.generate_min_items),
            ('maxItems', self.generate_max_items),
            ('uniqueItems', self.generate_unique_items),
            ('items', self.generate_items),
            ('minProperties', self.generate_min_properties),
            ('maxProperties', self.generate_max_properties),
            ('required', self.generate_required),
            # Check dependencies before properties generates default values.
            ('dependencies', self.generate_dependencies),
            ('properties', self.generate_properties),
            ('patternProperties', self.generate_pattern_properties),
            ('additionalProperties', self.generate_additional_properties),
        ))
        self._any_or_one_of_count = 0

    @property
    def global_state(self):
        res = super().global_state
        res['custom_formats'] = self._custom_formats
        return res

    def generate_type(self):
        """
        Validation of type. Can be one type or list of types.

        .. code-block:: python

            {'type': 'string'}
            {'type': ['string', 'number']}
        """
        types = enforce_list(self._definition['type'])
        try:
            python_types = ', '.join(JSON_TYPE_TO_PYTHON_TYPE[t] for t in types)
        except KeyError as exc:
            raise JsonSchemaDefinitionException('Unknown type: {}'.format(exc))

        extra = ''
        if ('number' in types or 'integer' in types) and 'boolean' not in types:
            extra = ' or isinstance({variable}, bool)'.format(variable=self._variable)

        with self.l('if not isinstance({variable}, ({})){}:', python_types, extra):
            self.exc('{name} must be {}', ' or '.join(types), rule='type')

    def generate_enum(self):
        """
        Means that only value specified in the enum is valid.

        .. code-block:: python

            {
                'enum': ['a', 'b'],
            }
        """
        enum = self._definition['enum']
        if not isinstance(enum, (list, tuple)):
            raise JsonSchemaDefinitionException('enum must be an array')
        with self.l('if {variable} not in {enum}:'):
            self.exc('{name} must be one of {}', self.e(enum), rule='enum')

    def generate_all_of(self):
        """
        Means that value have to be valid by all of those definitions. It's like put it in
        one big definition.

        .. code-block:: python

            {
                'allOf': [
                    {'type': 'number'},
                    {'minimum': 5},
                ],
            }

        Valid values for this definition are 5, 6, 7, ... but not 4 or 'abc' for example.
        """
        for definition_item in self._definition['allOf']:
            self.generate_func_code_block(definition_item, self._variable, self._variable_name, clear_variables=True)

    def generate_any_of(self):
        """
        Means that value have to be valid by any of those definitions. It can also be valid
        by all of them.

        .. code-block:: python

            {
                'anyOf': [
                    {'type': 'number', 'minimum': 10},
                    {'type': 'number', 'maximum': 5},
                ],
            }

        Valid values for this definition are 3, 4, 5, 10, 11, ... but not 8 for example.
        """
        self._any_or_one_of_count += 1
        count = self._any_or_one_of_count
        self.l('{variable}_any_of_count{count} = 0', count=count)
        for definition_item in self._definition['anyOf']:
            # When we know it's passing (at least once), we do not need to do another expensive try-except.
            with self.l('if not {variable}_any_of_count{count}:', count=count, optimize=False):
                with self.l('try:', optimize=False):
                    self.generate_func_code_block(definition_item, self._variable, self._variable_name, clear_variables=True)
                    self.l('{variable}_any_of_count{count} += 1', count=count)
                self.l('except JsonSchemaValueException: pass')

        with self.l('if not {variable}_any_of_count{count}:', count=count, optimize=False):
            self.exc('{name} cannot be validated by any definition', rule='anyOf')

    def generate_one_of(self):
        """
        Means that value have to be valid by only one of those definitions. It can't be valid
        by two or more of them.

        .. code-block:: python

            {
                'oneOf': [
                    {'type': 'number', 'multipleOf': 3},
                    {'type': 'number', 'multipleOf': 5},
                ],
            }

        Valid values for this definition are 3, 5, 6, ... but not 15 for example.
        """
        self._any_or_one_of_count += 1
        count = self._any_or_one_of_count
        self.l('{variable}_one_of_count{count} = 0', count=count)
        for definition_item in self._definition['oneOf']:
            # When we know it's failing (one of means exactly once), we do not need to do another expensive try-except.
            with self.l('if {variable}_one_of_count{count} < 2:', count=count, optimize=False):
                with self.l('try:', optimize=False):
                    self.generate_func_code_block(definition_item, self._variable, self._variable_name, clear_variables=True)
                    self.l('{variable}_one_of_count{count} += 1', count=count)
                self.l('except JsonSchemaValueException: pass')

        with self.l('if {variable}_one_of_count{count} != 1:', count=count):
            dynamic = '" (" + str({variable}_one_of_count{}) + " matches found)"'
            self.exc('{name} must be valid exactly by one definition', count, append_to_msg=dynamic, rule='oneOf')

    def generate_not(self):
        """
        Means that value have not to be valid by this definition.

        .. code-block:: python

            {'not': {'type': 'null'}}

        Valid values for this definition are 'hello', 42, {} ... but not None.

        Since draft 06 definition can be boolean. False means nothing, True
        means everything is invalid.
        """
        not_definition = self._definition['not']
        if not_definition is True:
            self.exc('{name} must not be there', rule='not')
        elif not_definition is False:
            return
        elif not not_definition:
            with self.l('if {}:', self._variable):
                self.exc('{name} must NOT match a disallowed definition', rule='not')
        else:
            with self.l('try:', optimize=False):
                self.generate_func_code_block(not_definition, self._variable, self._variable_name)
            self.l('except JsonSchemaValueException: pass')
            with self.l('else:'):
                self.exc('{name} must NOT match a disallowed definition', rule='not')

    def generate_min_length(self):
        with self.l('if isinstance({variable}, str):'):
            self.create_variable_with_length()
            if not isinstance(self._definition['minLength'], int):
                raise JsonSchemaDefinitionException('minLength must be a number')
            with self.l('if {variable}_len < {minLength}:'):
                self.exc('{name} must be longer than or equal to {minLength} characters', rule='minLength')

    def generate_max_length(self):
        with self.l('if isinstance({variable}, str):'):
            self.create_variable_with_length()
            if not isinstance(self._definition['maxLength'], int):
                raise JsonSchemaDefinitionException('maxLength must be a number')
            with self.l('if {variable}_len > {maxLength}:'):
                self.exc('{name} must be shorter than or equal to {maxLength} characters', rule='maxLength')

    def generate_pattern(self):
        with self.l('if isinstance({variable}, str):'):
            pattern = self._definition['pattern']
            safe_pattern = pattern.replace('\\', '\\\\').replace('"', '\\"')
            end_of_string_fixed_pattern = DOLLAR_FINDER.sub(r'\\Z', pattern)
            self._compile_regexps[pattern] = re.compile(end_of_string_fixed_pattern)
            with self.l('if not REGEX_PATTERNS[{}].search({variable}):', repr(pattern)):
                self.exc('{name} must match pattern {}', safe_pattern, rule='pattern')

    def generate_format(self):
        """
        Means that value have to be in specified format. For example date, email or other.

        .. code-block:: python

            {'format': 'email'}

        Valid value for this definition is user@example.com but not @username
        """
        if not self._use_formats:
            return
        with self.l('if isinstance({variable}, str):'):
            format_ = self._definition['format']
            # Checking custom formats - user is allowed to override default formats.
            if format_ in self._custom_formats:
                custom_format = self._custom_formats[format_]
                if isinstance(custom_format, str):
                    self._generate_format(format_, format_ + '_re_pattern', custom_format)
                else:
                    with self.l('if not custom_formats["{}"]({variable}):', format_):
                        self.exc('{name} must be {}', format_, rule='format')
            elif format_ in self.FORMAT_REGEXS:
                format_regex = self.FORMAT_REGEXS[format_]
                self._generate_format(format_, format_ + '_re_pattern', format_regex)
            # Format regex is used only in meta schemas.
            elif format_ == 'regex':
                with self.l('try:', optimize=False):
                    self.l('re.compile({variable})')
                with self.l('except Exception:'):
                    self.exc('{name} must be a valid regex', rule='format')
            else:
                raise JsonSchemaDefinitionException('Unknown format: {}'.format(format_))


    def _generate_format(self, format_name, regexp_name, regexp):
        if self._definition['format'] == format_name:
            if not regexp_name in self._compile_regexps:
                self._compile_regexps[regexp_name] = re.compile(regexp)
            with self.l('if not REGEX_PATTERNS["{}"].match({variable}):', regexp_name):
                self.exc('{name} must be {}', format_name, rule='format')

    def generate_minimum(self):
        with self.l('if isinstance({variable}, (int, float, Decimal)):'):
            if not isinstance(self._definition['minimum'], (int, float, decimal.Decimal)):
                raise JsonSchemaDefinitionException('minimum must be a number')
            if self._definition.get('exclusiveMinimum', False):
                with self.l('if {variable} <= {minimum}:'):
                    self.exc('{name} must be bigger than {minimum}', rule='minimum')
            else:
                with self.l('if {variable} < {minimum}:'):
                    self.exc('{name} must be bigger than or equal to {minimum}', rule='minimum')

    def generate_maximum(self):
        with self.l('if isinstance({variable}, (int, float, Decimal)):'):
            if not isinstance(self._definition['maximum'], (int, float, decimal.Decimal)):
                raise JsonSchemaDefinitionException('maximum must be a number')
            if self._definition.get('exclusiveMaximum', False):
                with self.l('if {variable} >= {maximum}:'):
                    self.exc('{name} must be smaller than {maximum}', rule='maximum')
            else:
                with self.l('if {variable} > {maximum}:'):
                    self.exc('{name} must be smaller than or equal to {maximum}', rule='maximum')

    def generate_multiple_of(self):
        with self.l('if isinstance({variable}, (int, float, Decimal)):'):
            if not isinstance(self._definition['multipleOf'], (int, float, decimal.Decimal)):
                raise JsonSchemaDefinitionException('multipleOf must be a number')
            # For proper multiplication check of floats we need to use decimals,
            # because for example 19.01 / 0.01 = 1901.0000000000002.
            if isinstance(self._definition['multipleOf'], float):
                self.l('quotient = Decimal(repr({variable})) / Decimal(repr({multipleOf}))')
            else:
                self.l('quotient = {variable} / {multipleOf}')
            with self.l('if int(quotient) != quotient:'):
                self.exc('{name} must be multiple of {multipleOf}', rule='multipleOf')

    def generate_min_items(self):
        self.create_variable_is_list()
        with self.l('if {variable}_is_list:'):
            if not isinstance(self._definition['minItems'], int):
                raise JsonSchemaDefinitionException('minItems must be a number')
            self.create_variable_with_length()
            with self.l('if {variable}_len < {minItems}:'):
                self.exc('{name} must contain at least {minItems} items', rule='minItems')

    def generate_max_items(self):
        self.create_variable_is_list()
        with self.l('if {variable}_is_list:'):
            if not isinstance(self._definition['maxItems'], int):
                raise JsonSchemaDefinitionException('maxItems must be a number')
            self.create_variable_with_length()
            with self.l('if {variable}_len > {maxItems}:'):
                self.exc('{name} must contain less than or equal to {maxItems} items', rule='maxItems')

    def generate_unique_items(self):
        """
        With Python 3.4 module ``timeit`` recommended this solutions:

        .. code-block:: python

            >>> timeit.timeit("len(x) > len(set(x))", "x=range(100)+range(100)", number=100000)
            0.5839540958404541
            >>> timeit.timeit("len({}.fromkeys(x)) == len(x)", "x=range(100)+range(100)", number=100000)
            0.7094449996948242
            >>> timeit.timeit("seen = set(); any(i in seen or seen.add(i) for i in x)", "x=range(100)+range(100)", number=100000)
            2.0819358825683594
            >>> timeit.timeit("np.unique(x).size == len(x)", "x=range(100)+range(100); import numpy as np", number=100000)
            2.1439831256866455
        """
        unique_definition = self._definition['uniqueItems']
        if not unique_definition:
            return

        self.create_variable_is_list()
        with self.l('if {variable}_is_list:'):
            self.l(
                'def fn(var): '
                'return frozenset(dict((k, fn(v)) '
                'for k, v in var.items()).items()) '
                'if hasattr(var, "items") else tuple(fn(v) '
                'for v in var) '
                'if isinstance(var, (dict, list)) else str(var) '
                'if isinstance(var, bool) else var')
            self.create_variable_with_length()
            with self.l('if {variable}_len > len(set(fn({variable}_x) for {variable}_x in {variable})):'):
                self.exc('{name} must contain unique items', rule='uniqueItems')

    def generate_items(self):
        """
        Means array is valid only when all items are valid by this definition.

        .. code-block:: python

            {
                'items': [
                    {'type': 'integer'},
                    {'type': 'string'},
                ],
            }

        Valid arrays are those with integers or strings, nothing else.

        Since draft 06 definition can be also boolean. True means nothing, False
        means everything is invalid.
        """
        items_definition = self._definition['items']
        if items_definition is True:
            return

        self.create_variable_is_list()
        with self.l('if {variable}_is_list:'):
            self.create_variable_with_length()
            if items_definition is False:
                with self.l('if {variable}:'):
                    self.exc('{name} must not be there', rule='items')
            elif isinstance(items_definition, list):
                for idx, item_definition in enumerate(items_definition):
                    with self.l('if {variable}_len > {}:', idx):
                        self.l('{variable}__{0} = {variable}[{0}]', idx)
                        self.generate_func_code_block(
                            item_definition,
                            '{}__{}'.format(self._variable, idx),
                            '{}[{}]'.format(self._variable_name, idx),
                        )
                    if self._use_default and isinstance(item_definition, dict) and 'default' in item_definition:
                        self.l('else: {variable}.append({})', repr(item_definition['default']))

                if 'additionalItems' in self._definition:
                    if self._definition['additionalItems'] is False:
                        with self.l('if {variable}_len > {}:', len(items_definition)):
                            self.exc('{name} must contain only specified items', rule='items')
                    else:
                        with self.l('for {variable}_x, {variable}_item in enumerate({variable}[{0}:], {0}):', len(items_definition)):
                            count = self.generate_func_code_block(
                                self._definition['additionalItems'],
                                '{}_item'.format(self._variable),
                                '{}[{{{}_x}}]'.format(self._variable_name, self._variable),
                            )
                            if count == 0:
                                self.l('pass')
            else:
                if items_definition:
                    with self.l('for {variable}_x, {variable}_item in enumerate({variable}):'):
                        count = self.generate_func_code_block(
                            items_definition,
                            '{}_item'.format(self._variable),
                            '{}[{{{}_x}}]'.format(self._variable_name, self._variable),
                        )
                        if count == 0:
                            self.l('pass')

    def generate_min_properties(self):
        self.create_variable_is_dict()
        with self.l('if {variable}_is_dict:'):
            if not isinstance(self._definition['minProperties'], int):
                raise JsonSchemaDefinitionException('minProperties must be a number')
            self.create_variable_with_length()
            with self.l('if {variable}_len < {minProperties}:'):
                self.exc('{name} must contain at least {minProperties} properties', rule='minProperties')

    def generate_max_properties(self):
        self.create_variable_is_dict()
        with self.l('if {variable}_is_dict:'):
            if not isinstance(self._definition['maxProperties'], int):
                raise JsonSchemaDefinitionException('maxProperties must be a number')
            self.create_variable_with_length()
            with self.l('if {variable}_len > {maxProperties}:'):
                self.exc('{name} must contain less than or equal to {maxProperties} properties', rule='maxProperties')

    def generate_required(self):
        self.create_variable_is_dict()
        with self.l('if {variable}_is_dict:'):
            if not isinstance(self._definition['required'], (list, tuple)):
                raise JsonSchemaDefinitionException('required must be an array')
            if len(self._definition['required']) != len(set(self._definition['required'])):
                raise JsonSchemaDefinitionException('required must contain unique elements')
            if not self._definition.get('additionalProperties', True):
                not_possible = [
                    prop
                    for prop in self._definition['required']
                    if
                        prop not in self._definition.get('properties', {})
                        and not any(re.search(regex, prop) for regex in self._definition.get('patternProperties', {}))
                ]
                if not_possible:
                    raise JsonSchemaDefinitionException('{}: items {} are required but not allowed'.format(self._variable, not_possible))
            self.l('{variable}__missing_keys = set({required}) - {variable}.keys()')
            with self.l('if {variable}__missing_keys:'):
                dynamic = 'str(sorted({variable}__missing_keys)) + " properties"'
                self.exc('{name} must contain ', self.e(self._definition['required']), rule='required', append_to_msg=dynamic)

    def generate_properties(self):
        """
        Means object with defined keys.

        .. code-block:: python

            {
                'properties': {
                    'key': {'type': 'number'},
                },
            }

        Valid object is containing key called 'key' and value any number.
        """
        self.create_variable_is_dict()
        with self.l('if {variable}_is_dict:'):
            self.create_variable_keys()
            for key, prop_definition in self._definition['properties'].items():
                key_name = re.sub(r'($[^a-zA-Z]|[^a-zA-Z0-9])', '', key)
                if not isinstance(prop_definition, (dict, bool)):
                    raise JsonSchemaDefinitionException('{}[{}] must be object'.format(self._variable, key_name))
                with self.l('if "{}" in {variable}_keys:', self.e(key)):
                    self.l('{variable}_keys.remove("{}")', self.e(key))
                    self.l('{variable}__{0} = {variable}["{1}"]', key_name, self.e(key))
                    self.generate_func_code_block(
                        prop_definition,
                        '{}__{}'.format(self._variable, key_name),
                        '{}.{}'.format(self._variable_name, self.e(key)),
                        clear_variables=True,
                    )
                if self._use_default and isinstance(prop_definition, dict) and 'default' in prop_definition:
                    self.l('else: {variable}["{}"] = {}', self.e(key), repr(prop_definition['default']))

    def generate_pattern_properties(self):
        """
        Means object with defined keys as patterns.

        .. code-block:: python

            {
                'patternProperties': {
                    '^x': {'type': 'number'},
                },
            }

        Valid object is containing key starting with a 'x' and value any number.
        """
        self.create_variable_is_dict()
        with self.l('if {variable}_is_dict:'):
            self.create_variable_keys()
            for pattern, definition in self._definition['patternProperties'].items():
                self._compile_regexps[pattern] = re.compile(pattern)
            with self.l('for {variable}_key, {variable}_val in {variable}.items():'):
                for pattern, definition in self._definition['patternProperties'].items():
                    with self.l('if REGEX_PATTERNS[{}].search({variable}_key):', repr(pattern)):
                        with self.l('if {variable}_key in {variable}_keys:'):
                            self.l('{variable}_keys.remove({variable}_key)')
                        self.generate_func_code_block(
                            definition,
                            '{}_val'.format(self._variable),
                            '{}.{{{}_key}}'.format(self._variable_name, self._variable),
                            clear_variables=True,
                        )

    def generate_additional_properties(self):
        """
        Means object with keys with values defined by definition.

        .. code-block:: python

            {
                'properties': {
                    'key': {'type': 'number'},
                }
                'additionalProperties': {'type': 'string'},
            }

        Valid object is containing key called 'key' and it's value any number and
        any other key with any string.
        """
        self.create_variable_is_dict()
        with self.l('if {variable}_is_dict:'):
            self.create_variable_keys()
            add_prop_definition = self._definition["additionalProperties"]
            if add_prop_definition is True or add_prop_definition == {}:
                return
            if add_prop_definition:
                properties_keys = list(self._definition.get("properties", {}).keys())
                with self.l('for {variable}_key in {variable}_keys:'):
                    with self.l('if {variable}_key not in {}:', properties_keys):
                        self.l('{variable}_value = {variable}.get({variable}_key)')
                        self.generate_func_code_block(
                            add_prop_definition,
                            '{}_value'.format(self._variable),
                            '{}.{{{}_key}}'.format(self._variable_name, self._variable),
                        )
            else:
                with self.l('if {variable}_keys:'):
                    self.exc('{name} must not contain "+str({variable}_keys)+" properties', rule='additionalProperties')

    def generate_dependencies(self):
        """
        Means when object has property, it needs to have also other property.

        .. code-block:: python

            {
                'dependencies': {
                    'bar': ['foo'],
                },
            }

        Valid object is containing only foo, both bar and foo or none of them, but not
        object with only bar.

        Since draft 06 definition can be boolean or empty array. True and empty array
        means nothing, False means that key cannot be there at all.
        """
        self.create_variable_is_dict()
        with self.l('if {variable}_is_dict:'):
            is_empty = True
            for key, values in self._definition["dependencies"].items():
                if values == [] or values is True:
                    continue
                is_empty = False
                with self.l('if "{}" in {variable}:', self.e(key)):
                    if values is False:
                        self.exc('{} in {name} must not be there', key, rule='dependencies')
                    elif isinstance(values, list):
                        for value in values:
                            with self.l('if "{}" not in {variable}:', self.e(value)):
                                self.exc('{name} missing dependency {} for {}', self.e(value), self.e(key), rule='dependencies')
                    else:
                        self.generate_func_code_block(values, self._variable, self._variable_name, clear_variables=True)
            if is_empty:
                self.l('pass')
