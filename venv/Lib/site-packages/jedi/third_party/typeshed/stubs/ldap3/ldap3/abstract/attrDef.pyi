from _typeshed import Incomplete

class AttrDef:
    name: Incomplete
    key: Incomplete
    validate: Incomplete
    pre_query: Incomplete
    post_query: Incomplete
    default: Incomplete
    dereference_dn: Incomplete
    description: Incomplete
    mandatory: Incomplete
    single_value: Incomplete
    oid_info: Incomplete
    other_names: Incomplete
    def __init__(
        self,
        name,
        key=None,
        validate=None,
        pre_query=None,
        post_query=None,
        default=...,
        dereference_dn=None,
        description=None,
        mandatory: bool = False,
        single_value=None,
        alias=None,
    ) -> None: ...
    def __eq__(self, other): ...
    def __lt__(self, other): ...
    def __hash__(self) -> int: ...
    def __setattr__(self, key: str, value) -> None: ...
