"""Mixin and helpers for participating in Resource Groups Tagging API operations."""

from __future__ import annotations

from collections.abc import Callable, Iterator
from dataclasses import dataclass, field
from typing import Any, ClassVar

from moto.core.base_backend import BackendDict
from moto.utilities.utils import get_partition


@dataclass
class TaggedResource:
    arn: str
    tags: dict[str, str]
    resource_type: str
    extra: dict[str, Any] = field(default_factory=dict)


class TaggableResourcesMixin:
    """Mix into a service backend to participate in Resource Groups Tagging API.

    Subclasses must set `SERVICE_NAMESPACE` (the AWS service identifier used in
    ARNs, e.g. "rds") and override `iter_tagged_resources`. `tag_resource`
    and `untag_resource` only need overriding if the service supports the
    cross-service TagResources / UntagResources operations.
    """

    SERVICE_NAMESPACE: ClassVar[str] = ""

    def __init_subclass__(cls, **kwargs: Any) -> None:
        super().__init_subclass__(**kwargs)
        if not getattr(cls, "SERVICE_NAMESPACE", ""):
            raise TypeError(
                f"{cls.__name__} mixes in TaggableResourcesMixin but does not set SERVICE_NAMESPACE"
            )

    def iter_tagged_resources(self) -> Iterator[TaggedResource]:
        """Yield every taggable resource owned by this backend.

        The central code applies tag filters and resource-type filters, so implementations
        should yield resources unconditionally (even those with no tags).
        """
        return iter(())

    def owns_arn(self, arn: str) -> bool:
        """Return True if this backend should handle tag/untag for `arn`.

        Default matches `arn:<partition>:<SERVICE_NAME>:`. Override when the ARN scheme
        doesn't match the service name (e.g. ELBv2 uses `elasticloadbalancing`).
        """
        return arn.startswith(f"arn:{self.partition}:{self.SERVICE_NAMESPACE}:")  # type: ignore[attr-defined]

    def tag_resource(self, arn: str, tags: dict[str, str]) -> None:
        raise NotImplementedError(
            f"{type(self).__name__} does not implement tag_resource"
        )

    def untag_resource(self, arn: str, tag_keys: list[str]) -> None:
        raise NotImplementedError(
            f"{type(self).__name__} does not implement untag_resource"
        )


def iter_taggable_backends(
    account_id: str, region_name: str
) -> Iterator[TaggableResourcesMixin]:
    """Yield every in-use TaggableResourcesMixin backend for the given account/region.

    Walks a list of BackendDicts that have at least one account-specific backend created.
    Services the user has never touched are skipped (as they have no resources to report).
    """
    for backend_dict in list(BackendDict._instances):  # type: ignore[misc]
        if not isinstance(backend_dict.backend, type):
            continue
        if not issubclass(backend_dict.backend, TaggableResourcesMixin):
            continue
        if account_id not in backend_dict:
            continue
        account_specific_backend = backend_dict[account_id]
        if (
            region_name in account_specific_backend
            or region_name in account_specific_backend.regions
        ):
            yield account_specific_backend[region_name]
            continue
        # Yield global services keyed by partition (e.g. S3).
        partition = get_partition(region_name)
        if partition in account_specific_backend.keys():
            yield account_specific_backend[partition]


def match_resource_type(
    resource_type: str, resource_type_filters: list[str] | None
) -> bool:
    """Return True if `resource_type` (e.g. "rds:cluster") satisfies the filters.

    A filter of "rds" matches every "rds:*" resource type; a filter of
    "rds:cluster" matches only that exact resource type. If no filters
    are provided, all resource types match.
    """
    if not resource_type_filters:
        return True
    service = resource_type.split(":", 1)[0]
    for f in resource_type_filters:
        if f == resource_type or f == service:
            return True
    return False


def make_tag_matcher(
    tag_filters: list[dict[str, Any]] | None,
) -> Callable[[dict[str, str]], bool]:
    """Return a predicate matching the GetResources TagFilters semantics.

    A filter is satisfied when:
      - if no Values are given: any tag with a matching Key is present
      - if Values are given: a tag with a matching Key and a Value in Values is present.
    All filters must be satisfied (AND when multiple filters, OR when multiple values).
    """
    filters = list(tag_filters or [])

    def matches(tags: dict[str, str]) -> bool:
        if not filters:
            return True
        tags_by_key: dict[str, list[str]] = {}
        for tag_key, tag_value in tags.items():
            tags_by_key.setdefault(tag_key, []).append(tag_value)
        for f in filters:
            key = f["Key"]
            values = f.get("Values") or []
            if key not in tags_by_key:
                return False
            if values and not any(v in tags_by_key[key] for v in values):
                return False
        return True

    return matches
