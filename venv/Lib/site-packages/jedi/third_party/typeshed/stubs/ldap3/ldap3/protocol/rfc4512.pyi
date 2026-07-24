from _typeshed import Incomplete

def constant_to_class_kind(value): ...
def constant_to_attribute_usage(value): ...
def attribute_usage_to_constant(value): ...
def quoted_string_to_list(quoted_string): ...
def oids_string_to_list(oid_string): ...
def extension_to_tuple(extension_string): ...
def list_to_string(list_object): ...

class BaseServerInfo:
    raw: Incomplete
    def __init__(self, raw_attributes) -> None: ...
    @classmethod
    def from_json(cls, json_definition, schema=None, custom_formatter=None): ...
    @classmethod
    def from_file(cls, target, schema=None, custom_formatter=None): ...
    def to_file(self, target, indent: int = 4, sort: bool = True) -> None: ...
    def to_json(self, indent: int = 4, sort: bool = True): ...

class DsaInfo(BaseServerInfo):
    alt_servers: Incomplete
    naming_contexts: Incomplete
    supported_controls: Incomplete
    supported_extensions: Incomplete
    supported_features: Incomplete
    supported_ldap_versions: Incomplete
    supported_sasl_mechanisms: Incomplete
    vendor_name: Incomplete
    vendor_version: Incomplete
    schema_entry: Incomplete
    other: Incomplete
    def __init__(self, attributes, raw_attributes) -> None: ...

class SchemaInfo(BaseServerInfo):
    schema_entry: Incomplete
    create_time_stamp: Incomplete
    modify_time_stamp: Incomplete
    attribute_types: Incomplete
    object_classes: Incomplete
    matching_rules: Incomplete
    matching_rule_uses: Incomplete
    dit_content_rules: Incomplete
    dit_structure_rules: Incomplete
    name_forms: Incomplete
    ldap_syntaxes: Incomplete
    other: Incomplete
    def __init__(self, schema_entry, attributes, raw_attributes) -> None: ...
    def is_valid(self): ...

class BaseObjectInfo:
    oid: Incomplete
    name: Incomplete
    description: Incomplete
    obsolete: Incomplete
    extensions: Incomplete
    experimental: Incomplete
    raw_definition: Incomplete
    def __init__(
        self, oid=None, name=None, description=None, obsolete: bool = False, extensions=None, experimental=None, definition=None
    ) -> None: ...
    @property
    def oid_info(self): ...
    @classmethod
    def from_definition(cls, definitions): ...

class MatchingRuleInfo(BaseObjectInfo):
    syntax: Incomplete
    def __init__(
        self,
        oid=None,
        name=None,
        description=None,
        obsolete: bool = False,
        syntax=None,
        extensions=None,
        experimental=None,
        definition=None,
    ) -> None: ...

class MatchingRuleUseInfo(BaseObjectInfo):
    apply_to: Incomplete
    def __init__(
        self,
        oid=None,
        name=None,
        description=None,
        obsolete: bool = False,
        apply_to=None,
        extensions=None,
        experimental=None,
        definition=None,
    ) -> None: ...

class ObjectClassInfo(BaseObjectInfo):
    superior: Incomplete
    kind: Incomplete
    must_contain: Incomplete
    may_contain: Incomplete
    def __init__(
        self,
        oid=None,
        name=None,
        description=None,
        obsolete: bool = False,
        superior=None,
        kind=None,
        must_contain=None,
        may_contain=None,
        extensions=None,
        experimental=None,
        definition=None,
    ) -> None: ...

class AttributeTypeInfo(BaseObjectInfo):
    superior: Incomplete
    equality: Incomplete
    ordering: Incomplete
    substring: Incomplete
    syntax: Incomplete
    min_length: Incomplete
    single_value: Incomplete
    collective: Incomplete
    no_user_modification: Incomplete
    usage: Incomplete
    mandatory_in: Incomplete
    optional_in: Incomplete
    def __init__(
        self,
        oid=None,
        name=None,
        description=None,
        obsolete: bool = False,
        superior=None,
        equality=None,
        ordering=None,
        substring=None,
        syntax=None,
        min_length=None,
        single_value: bool = False,
        collective: bool = False,
        no_user_modification: bool = False,
        usage=None,
        extensions=None,
        experimental=None,
        definition=None,
    ) -> None: ...

class LdapSyntaxInfo(BaseObjectInfo):
    def __init__(self, oid=None, description=None, extensions=None, experimental=None, definition=None) -> None: ...

class DitContentRuleInfo(BaseObjectInfo):
    auxiliary_classes: Incomplete
    must_contain: Incomplete
    may_contain: Incomplete
    not_contains: Incomplete
    def __init__(
        self,
        oid=None,
        name=None,
        description=None,
        obsolete: bool = False,
        auxiliary_classes=None,
        must_contain=None,
        may_contain=None,
        not_contains=None,
        extensions=None,
        experimental=None,
        definition=None,
    ) -> None: ...

class DitStructureRuleInfo(BaseObjectInfo):
    superior: Incomplete
    name_form: Incomplete
    def __init__(
        self,
        oid=None,
        name=None,
        description=None,
        obsolete: bool = False,
        name_form=None,
        superior=None,
        extensions=None,
        experimental=None,
        definition=None,
    ) -> None: ...

class NameFormInfo(BaseObjectInfo):
    object_class: Incomplete
    must_contain: Incomplete
    may_contain: Incomplete
    def __init__(
        self,
        oid=None,
        name=None,
        description=None,
        obsolete: bool = False,
        object_class=None,
        must_contain=None,
        may_contain=None,
        extensions=None,
        experimental=None,
        definition=None,
    ) -> None: ...
