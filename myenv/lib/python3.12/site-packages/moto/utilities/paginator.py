import inspect
from copy import deepcopy
from functools import wraps
from typing import TYPE_CHECKING, Any, Callable, Dict, List, Optional, Tuple, TypeVar

from botocore.paginate import TokenDecoder, TokenEncoder

from moto.core.exceptions import InvalidToken

if TYPE_CHECKING:
    from typing_extensions import ParamSpec, Protocol

    P1 = ParamSpec("P1")
    P2 = ParamSpec("P2")
else:
    Protocol = object

T = TypeVar("T")


class GenericFunction(Protocol):
    def __call__(
        self, func: "Callable[P1, List[T]]"
    ) -> "Callable[P2, Tuple[List[T], Optional[str]]]": ...


def paginate(pagination_model: Dict[str, Any]) -> GenericFunction:
    def pagination_decorator(
        func: Callable[..., List[T]],
    ) -> Callable[..., Tuple[List[T], Optional[str]]]:
        @wraps(func)
        def pagination_wrapper(*args: Any, **kwargs: Any) -> Any:  # type: ignore
            method = func.__name__
            model = pagination_model
            pagination_config = model.get(method)
            if not pagination_config:
                raise ValueError(f"No pagination config for backend method: {method}")
            # Get the pagination arguments, to be used by the paginator
            next_token_name = pagination_config.get("input_token", "next_token")
            limit_name = pagination_config.get("limit_key")
            input_token = kwargs.get(next_token_name)
            limit = kwargs.get(limit_name, None)
            # Remove pagination arguments from our input kwargs
            # We need this to verify that our input kwargs are the same across invocations
            #   list_all(service="x")                   next_token = "a"
            #   list_all(service="x", next_token="a") ==> Works fine
            #   list_all(service="y", next_token="a") ==> Should throw an error, as the input_kwargs are different
            input_kwargs = deepcopy(kwargs)
            input_kwargs.pop(next_token_name, None)
            input_kwargs.pop(limit_name, None)
            fail_on_invalid_token = pagination_config.get("fail_on_invalid_token", True)
            paginator = Paginator(
                max_results=limit,
                max_results_default=pagination_config.get("limit_default"),
                starting_token=input_token,
                unique_attribute=pagination_config.get("unique_attribute"),
                param_values_to_check=input_kwargs,
                fail_on_invalid_token=fail_on_invalid_token,
            )

            # Determine which parameters to pass
            (arg_names, _, has_kwargs, _, _, _, _) = inspect.getfullargspec(func)
            # If the target-func expects `**kwargs`, we can pass everything
            if not has_kwargs:
                # If the target-function does not expect the next_token/limit, do not pass it
                if next_token_name not in arg_names:
                    kwargs.pop(next_token_name, None)
                if limit_name not in arg_names:
                    kwargs.pop(limit_name, None)

            results = func(*args, **kwargs)
            return paginator.paginate(results)

        return pagination_wrapper

    return pagination_decorator


class Paginator:
    def __init__(
        self,
        max_results: Any = None,
        max_results_default: Any = None,
        starting_token: Any = None,
        unique_attribute: Any = None,
        param_values_to_check: Any = None,
        fail_on_invalid_token: bool = True,
    ):
        self._max_results = max_results if max_results else max_results_default
        self._starting_token = starting_token
        self._unique_attributes = (
            unique_attribute
            if isinstance(unique_attribute, list)
            else [unique_attribute]
        )
        self._param_values_to_check = param_values_to_check
        self._fail_on_invalid_token = fail_on_invalid_token
        self._token_encoder = TokenEncoder()
        self._token_decoder = TokenDecoder()
        self._param_checksum = self._calculate_parameter_checksum()
        self._parsed_token = self._parse_starting_token()

    def _parse_starting_token(self) -> Optional[Dict[str, Any]]:
        if self._starting_token is None:
            return None
        # The starting token is a dict passed as a base64 encoded string.
        next_token = self._starting_token
        try:
            next_token = self._token_decoder.decode(next_token)
        except (ValueError, TypeError, UnicodeDecodeError):
            self._raise_exception_if_required(next_token)
            return None
        if next_token.get("parameterChecksum") != self._param_checksum:
            raise InvalidToken(f"Input inconsistent with page token: {str(next_token)}")
        return next_token

    def _raise_exception_if_required(self, token: Optional[str]) -> None:
        if self._fail_on_invalid_token:
            if isinstance(self._fail_on_invalid_token, type):
                # we need to raise a custom exception
                func_info = inspect.getfullargspec(self._fail_on_invalid_token)
                arg_names, _, _, _, _, _, _ = func_info
                # arg_names == [self] or [self, token_argument_that_can_have_any_name]
                requires_token_arg = len(arg_names) > 1
                if requires_token_arg:
                    raise self._fail_on_invalid_token(token)
                else:
                    raise self._fail_on_invalid_token()
            raise InvalidToken("Invalid token")

    def _calculate_parameter_checksum(self) -> int:
        def freeze(o: Any) -> Any:
            if not o:
                return None
            if isinstance(o, dict):
                return frozenset({k: freeze(v) for k, v in o.items()}.items())

            if isinstance(o, (list, tuple, set)):
                return tuple([freeze(v) for v in o])

            return o

        return hash(freeze(self._param_values_to_check))

    def _check_predicate(self, item: Any) -> bool:
        if self._parsed_token is None:
            return False
        unique_attributes = self._parsed_token["uniqueAttributes"]
        predicate_values = unique_attributes.split("|")
        for index, attr in enumerate(self._unique_attributes):
            curr_val = (
                item[attr] if isinstance(item, dict) else getattr(item, attr, None)
            )
            if not str(curr_val) == predicate_values[index]:
                return False
        return True

    def _build_next_token(self, next_item: Any) -> str:
        token_dict: Dict[str, Any] = {}
        if self._param_checksum:
            token_dict["parameterChecksum"] = self._param_checksum
        range_keys = []
        for attr in self._unique_attributes:
            if isinstance(next_item, dict):
                range_keys.append(str(next_item[attr]))
            else:
                range_keys.append(str(getattr(next_item, attr)))
        token_dict["uniqueAttributes"] = "|".join(range_keys)
        return self._token_encoder.encode(token_dict)

    def paginate(self, results: List[Any]) -> Tuple[List[Any], Optional[str]]:
        index_start = 0
        if self._starting_token:
            try:
                index_start = next(
                    index
                    for (index, result) in enumerate(results)
                    if self._check_predicate(result)
                )
            except StopIteration:
                if self._fail_on_invalid_token:
                    raise InvalidToken("Resource not found!")
                else:
                    return [], None

        index_end = index_start + self._max_results
        if index_end > len(results):
            index_end = len(results)

        results_page = results[index_start:index_end]

        next_token = None
        if results_page and index_end < len(results):
            last_resource_on_this_page = results[index_end]
            next_token = self._build_next_token(last_resource_on_this_page)
        return results_page, next_token
