import abc
import logging
import threading
from typing import Any, Callable, Dict, List, TypedDict, Union

from moto.moto_api._internal import mock_random

log = logging.getLogger(__name__)

ExistingIds = Union[List[str], None]
Tags = Union[Dict[str, str], List[Dict[str, str]], None]

# Custom resource tag to override the generated resource ID.
TAG_KEY_CUSTOM_ID = "_custom_id_"


class IdSourceContext(TypedDict, total=False):
    resource_identifier: "ResourceIdentifier"
    tags: Tags
    existing_ids: ExistingIds


class ResourceIdentifier(abc.ABC):
    """
    Base class for resource identifiers. When implementing a new resource, it is important to set
    the service and resource as they will be used to create the unique identifier for that resource.

    It is recommended to implement the `generate` method using functions decorated with `@moto_id`.
    This will ensure that your resource can be assigned a custom id.
    """

    service: str
    resource: str

    def __init__(self, account_id: str, region: str, name: str):
        self.account_id = account_id
        self.region = region
        self.name = name or ""

    @abc.abstractmethod
    def generate(self, existing_ids: ExistingIds = None, tags: Tags = None) -> str:
        """Method to generate a resource id"""

    @property
    def unique_identifier(self) -> str:
        return ".".join(
            [self.account_id, self.region, self.service, self.resource, self.name]
        )

    def __str__(self) -> str:
        return self.unique_identifier


class MotoIdManager:
    """class to manage custom ids. Do not create instance and instead
    use the `id_manager` instance created below."""

    _custom_ids: Dict[str, str]
    _id_sources: List[Callable[[IdSourceContext], Union[str, None]]]

    _lock: threading.RLock

    def __init__(self) -> None:
        self._custom_ids = {}
        self._lock = threading.RLock()
        self._id_sources = []

        self.add_id_source(self.get_id_from_tags)
        self.add_id_source(self.get_custom_id_from_context)

    def get_custom_id(
        self, resource_identifier: ResourceIdentifier
    ) -> Union[str, None]:
        # retrieves a custom_id for a resource. Returns None if no id were registered
        # that matches the `resource_identifier`
        return self._custom_ids.get(resource_identifier.unique_identifier)

    def set_custom_id(
        self, resource_identifier: ResourceIdentifier, custom_id: str
    ) -> None:
        # Do not set a custom_id for a resource no value was found for the name
        if not resource_identifier.name:
            return
        with self._lock:
            self._custom_ids[resource_identifier.unique_identifier] = custom_id

    def unset_custom_id(self, resource_identifier: ResourceIdentifier) -> None:
        # removes a set custom_id for a resource
        with self._lock:
            self._custom_ids.pop(resource_identifier.unique_identifier, None)

    def add_id_source(
        self, id_source: Callable[[IdSourceContext], Union[str, None]]
    ) -> None:
        self._id_sources.append(id_source)

    @staticmethod
    def get_id_from_tags(id_source_context: IdSourceContext) -> Union[str, None]:
        if not (tags := id_source_context.get("tags")):
            return None

        if isinstance(tags, dict):
            return tags.get(TAG_KEY_CUSTOM_ID)

        if isinstance(tags, list):
            return next(
                (
                    tag.get("Value")
                    for tag in tags
                    if tag.get("Key") == TAG_KEY_CUSTOM_ID
                ),
                None,
            )

    def get_custom_id_from_context(
        self, id_source_context: IdSourceContext
    ) -> Union[str, None]:
        # retrieves a custom_id for a resource. Returns None
        if resource_identifier := id_source_context.get("resource_identifier"):
            return self.get_custom_id(resource_identifier)
        return None

    def find_id_from_sources(
        self, id_source_context: IdSourceContext
    ) -> Union[str, None]:
        existing_ids = id_source_context.get("existing_ids") or []
        for id_source in self._id_sources:
            if found_id := id_source(id_source_context):
                if found_id in existing_ids:
                    log.debug(
                        f"Found id {found_id} for resource {id_source_context.get('resource_identifier')}, "
                        "but a resource already exists with this id."
                    )
                else:
                    return found_id

        return None


moto_id_manager = MotoIdManager()


def moto_id(fn: Callable[..., str]) -> Callable[..., str]:
    """
    Decorator for helping in creation of static ids.

    The decorated function should accept the following parameters

    :param resource_identifier
    :param existing_ids
        If provided, we will omit returning a custom id if it is already on the list
    :param tags
        If provided will look for a tag named `_custom_id_`. This will take precedence over registered custom ids
    """

    def _wrapper(
        resource_identifier: ResourceIdentifier,
        existing_ids: ExistingIds = None,
        tags: Tags = None,
        **kwargs: Dict[str, Any],
    ) -> str:
        if resource_identifier and (
            found_id := moto_id_manager.find_id_from_sources(
                IdSourceContext(
                    resource_identifier=resource_identifier,
                    existing_ids=existing_ids,
                    tags=tags,
                )
            )
        ):
            return found_id

        return fn(
            resource_identifier=resource_identifier,
            existing_ids=existing_ids,
            tags=tags,
            **kwargs,
        )

    return _wrapper


@moto_id
def generate_str_id(  # type: ignore
    resource_identifier: ResourceIdentifier,
    existing_ids: ExistingIds = None,
    tags: Tags = None,
    length: int = 20,
    include_digits: bool = True,
    lower_case: bool = False,
) -> str:
    return mock_random.get_random_string(length, include_digits, lower_case)
