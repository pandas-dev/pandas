from collections.abc import Iterator
from typing import Any

from moto.core.base_backend import BackendDict, BaseBackend
from moto.core.resource_tagging import (
    TaggableResourcesMixin,
    TaggedResource,
    iter_taggable_backends,
    make_tag_matcher,
    match_resource_type,
)
from moto.moto_api._internal import mock_random
from moto.resourcegroupstaggingapi.exceptions import (
    ResourceGroupsTaggingAPIError as RESTError,
)


class ResourceGroupsTaggingAPIBackend(BaseBackend):
    def __init__(self, region_name: str, account_id: str):
        super().__init__(region_name, account_id)

        self._pages: dict[str, Any] = {}
        # Like 'someuuid': {'gen': <generator>, 'misc': None}
        # Misc is there for peeking from a generator and it cant
        # fit in the current request. As we only store generators
        # there is really no point cleaning up

    def _get_backend_for_resource(
        self, resource_arn: str
    ) -> TaggableResourcesMixin | None:
        backend = next(
            (
                b
                for b in iter_taggable_backends(self.account_id, self.region_name)
                if b.owns_arn(resource_arn)
            ),
            None,
        )
        return backend

    def _get_resources_generator(
        self,
        tag_filters: list[dict[str, Any]] | None = None,
        resource_type_filters: list[str] | None = None,
        resource_arn_list: list[str] | None = None,
    ) -> Iterator[TaggedResource]:
        tag_matcher = make_tag_matcher(tag_filters)
        for backend in iter_taggable_backends(self.account_id, self.region_name):
            for resource in backend.iter_tagged_resources():
                # https://docs.aws.amazon.com/resourcegroupstagging/latest/APIReference/API_GetResources.html
                # According to the docs, GetResources returns all "tagged or previously tagged" resources.
                # This means that even though generally resources without tags are not returned, a resource
                # that had been tagged and then had all of its tags deleted would be returned.
                # Historically, Moto has not supported this behavior and has been inconsistent across
                # service backends with respect to returning resources without tags (some do, some don't).
                # The extra "include_untagged" here allows backends to decide whether resources without
                # tags should be returned.
                # TODO: apply consistent behavior for untagged resources that matches real AWS.
                if not resource.tags and not resource.extra.get("include_untagged"):
                    continue
                if not match_resource_type(
                    resource.resource_type, resource_type_filters
                ):
                    continue
                if not tag_matcher(resource.tags):
                    continue
                if resource_arn_list and resource.arn not in resource_arn_list:
                    continue
                yield resource

    def _get_tag_keys_generator(self) -> Iterator[str]:
        for backend in iter_taggable_backends(self.account_id, self.region_name):
            for resource in backend.iter_tagged_resources():
                yield from resource.tags.keys()

    def _get_tag_values_generator(self, tag_key: str) -> Iterator[str]:
        for backend in iter_taggable_backends(self.account_id, self.region_name):
            for resource in backend.iter_tagged_resources():
                for key, value in resource.tags.items():
                    if key == tag_key:
                        yield value

    def get_resources(
        self,
        pagination_token: str | None = None,
        resources_per_page: int = 50,
        tags_per_page: int = 100,
        tag_filters: list[dict[str, Any]] | None = None,
        resource_type_filters: list[str] | None = None,
        resource_arn_list: list[str] | None = None,
    ) -> tuple[str | None, list[TaggedResource]]:
        # If we have a token, go and find the respective generator, or error
        if pagination_token:
            if pagination_token not in self._pages:
                raise RESTError(
                    "PaginationTokenExpiredException", "Token does not exist"
                )

            generator = self._pages[pagination_token]["gen"]
            left_over = self._pages[pagination_token]["misc"]
        else:
            generator = self._get_resources_generator(
                tag_filters=tag_filters,
                resource_type_filters=resource_type_filters,
                resource_arn_list=resource_arn_list,
            )
            left_over = None

        result = []
        current_tags = 0
        current_resources = 0
        if left_over:
            result.append(left_over)
            current_resources += 1

        try:
            while True:
                next_item = next(generator)
                resource_tags = len(next_item.tags)

                if current_resources >= resources_per_page:
                    break
                if current_tags + resource_tags >= tags_per_page:
                    break

                current_resources += 1
                current_tags += resource_tags

                result.append(next_item)

        except StopIteration:
            # Finished generator before invalidating page limiting constraints
            return None, result

        # Didn't hit StopIteration so there's stuff left in generator
        new_token = str(mock_random.uuid4())
        self._pages[new_token] = {"gen": generator, "misc": next_item}

        # Token used up, might as well bin now, if you call it again you're an idiot
        if pagination_token:
            del self._pages[pagination_token]

        return new_token, result

    def get_tag_keys(
        self, pagination_token: str | None = None
    ) -> tuple[str | None, list[str]]:
        if pagination_token:
            if pagination_token not in self._pages:
                raise RESTError(
                    "PaginationTokenExpiredException", "Token does not exist"
                )

            generator = self._pages[pagination_token]["gen"]
            left_over = self._pages[pagination_token]["misc"]
        else:
            generator = self._get_tag_keys_generator()
            left_over = None

        result = []
        current_tags = 0
        if left_over:
            result.append(left_over)
            current_tags += 1

        try:
            while True:
                # Generator format: ['tag', 'tag', 'tag', ...]
                next_item = next(generator)

                if current_tags + 1 >= 128:
                    break

                current_tags += 1

                result.append(next_item)

        except StopIteration:
            # Finished generator before invalidating page limiting constraints
            return None, result

        # Didn't hit StopIteration so there's stuff left in generator
        new_token = str(mock_random.uuid4())
        self._pages[new_token] = {"gen": generator, "misc": next_item}

        # Token used up, might as well bin now, if you call it again your an idiot
        if pagination_token:
            del self._pages[pagination_token]

        return new_token, result

    def get_tag_values(
        self, pagination_token: str | None, key: str
    ) -> tuple[str | None, list[str]]:
        if pagination_token:
            if pagination_token not in self._pages:
                raise RESTError(
                    "PaginationTokenExpiredException", "Token does not exist"
                )

            generator = self._pages[pagination_token]["gen"]
            left_over = self._pages[pagination_token]["misc"]
        else:
            generator = self._get_tag_values_generator(key)
            left_over = None

        result = []
        current_tags = 0
        if left_over:
            result.append(left_over)
            current_tags += 1

        try:
            while True:
                # Generator format: ['value', 'value', 'value', ...]
                next_item = next(generator)

                if current_tags + 1 >= 128:
                    break

                current_tags += 1

                result.append(next_item)

        except StopIteration:
            # Finished generator before invalidating page limiting constraints
            return None, result

        # Didn't hit StopIteration so there's stuff left in generator
        new_token = str(mock_random.uuid4())
        self._pages[new_token] = {"gen": generator, "misc": next_item}

        # Token used up, might as well bin now, if you call it again your an idiot
        if pagination_token:
            del self._pages[pagination_token]

        return new_token, result

    def tag_resources(
        self, resource_arns: list[str], tags: dict[str, str]
    ) -> dict[str, dict[str, Any]]:
        missing_resources = []
        missing_error: dict[str, Any] = {
            "StatusCode": 404,
            "ErrorCode": "InternalServiceException",
            "ErrorMessage": "Service not yet supported",
        }

        for arn in resource_arns:
            backend_for_arn = self._get_backend_for_resource(arn)
            if backend_for_arn is None:
                missing_resources.append(arn)
                continue
            try:
                backend_for_arn.tag_resource(arn, tags)
            except NotImplementedError:
                missing_resources.append(arn)
        return dict.fromkeys(missing_resources, missing_error)

    def untag_resources(
        self, resource_arn_list: list[str], tag_keys: list[str]
    ) -> dict[str, dict[str, Any]]:
        missing_resources = []
        missing_error: dict[str, Any] = {
            "StatusCode": 404,
            "ErrorCode": "InternalServiceException",
            "ErrorMessage": "Service not yet supported",
        }

        for arn in resource_arn_list:
            backend_for_arn = self._get_backend_for_resource(arn)
            if backend_for_arn is None:
                missing_resources.append(arn)
                continue
            try:
                backend_for_arn.untag_resource(arn, tag_keys)
            except NotImplementedError:
                missing_resources.append(arn)

        return dict.fromkeys(missing_resources, missing_error)


resourcegroupstaggingapi_backends = BackendDict(
    ResourceGroupsTaggingAPIBackend, "resourcegroupstaggingapi"
)
