from _typeshed import Incomplete
from collections.abc import Iterable
from typing import TypedDict, type_check_only

from braintree.configuration import Configuration
from braintree.environment import Environment
from braintree.util.http import Http

@type_check_only
class _Extension(TypedDict):
    errorClass: Incomplete
    legacyCode: int | None

@type_check_only
class _Error(TypedDict):
    attribute: str | None
    code: int | None
    message: str | None
    extensions: _Extension | None

@type_check_only
class _ValidationErrors(TypedDict):
    errors: Iterable[_Error]

@type_check_only
class _Response(TypedDict):
    errors: Iterable[_Error] | None

class GraphQLClient(Http):
    @staticmethod
    def raise_exception_for_graphql_error(response: _Response) -> None: ...
    graphql_headers: dict[str, str]
    def __init__(self, config: Configuration | None = None, environment: Environment | None = None) -> None: ...
    def query(self, definition, variables=None, operation_name=None): ...
    @staticmethod
    def get_validation_errors(response) -> _ValidationErrors | None: ...
    @staticmethod
    def get_validation_error_code(error: _Error) -> int | None: ...
