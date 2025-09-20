import re
from collections import defaultdict
from typing import Dict, List, Optional

from ..exceptions import (
    InvalidParameterValueErrorTagNull,
    TagLimitExceeded,
)
from ..utils import (
    EC2_PREFIX_TO_RESOURCE,
    get_prefix,
    simple_aws_filter_to_re,
)


class TagBackend:
    VALID_TAG_FILTERS = ["key", "resource-id", "resource-type", "value"]

    def __init__(self) -> None:
        self.tags: Dict[str, Dict[str, str]] = defaultdict(dict)

    def create_tags(self, resource_ids: List[str], tags: Dict[str, str]) -> bool:
        if None in set([tags[tag] for tag in tags]):
            raise InvalidParameterValueErrorTagNull()
        for resource_id in resource_ids:
            if resource_id in self.tags:
                if (
                    len(self.tags[resource_id])
                    + len(
                        [
                            tag
                            for tag in tags
                            if not tag.startswith("aws:")
                            and tag not in self.tags[resource_id]
                        ]
                    )
                    > 50
                ):
                    raise TagLimitExceeded()
            elif len([tag for tag in tags if not tag.startswith("aws:")]) > 50:
                raise TagLimitExceeded()
        for resource_id in resource_ids:
            for tag in tags:
                self.tags[resource_id][tag] = tags[tag]
        return True

    def delete_tags(self, resource_ids: List[str], tags: Dict[str, str]) -> bool:
        for resource_id in resource_ids:
            for tag in tags:
                if tag in self.tags[resource_id]:
                    if tags[tag] is None:
                        self.tags[resource_id].pop(tag)
                    elif tags[tag] == self.tags[resource_id][tag]:
                        self.tags[resource_id].pop(tag)
        return True

    def describe_tags(
        self, filters: Optional[Dict[str, List[str]]] = None
    ) -> List[Dict[str, str]]:
        results = []
        key_filters = []
        resource_id_filters = []
        resource_type_filters = []
        value_filters = []
        if filters is not None:
            for tag_filter in filters:
                if tag_filter in self.VALID_TAG_FILTERS:
                    if tag_filter == "key":
                        for value in filters[tag_filter]:
                            key_filters.append(
                                re.compile(simple_aws_filter_to_re(value))
                            )
                    if tag_filter == "resource-id":
                        for value in filters[tag_filter]:
                            resource_id_filters.append(
                                re.compile(simple_aws_filter_to_re(value))
                            )
                    if tag_filter == "resource-type":
                        for value in filters[tag_filter]:
                            resource_type_filters.append(value)
                    if tag_filter == "value":
                        for value in filters[tag_filter]:
                            value_filters.append(
                                re.compile(simple_aws_filter_to_re(value))
                            )
        for resource_id, tags in self.tags.copy().items():
            for key, value in tags.items():
                add_result = False
                if filters is None:
                    add_result = True
                else:
                    key_pass = False
                    id_pass = False
                    type_pass = False
                    value_pass = False
                    if key_filters:
                        for pattern in key_filters:
                            if pattern.match(key) is not None:
                                key_pass = True
                    else:
                        key_pass = True
                    if resource_id_filters:
                        for pattern in resource_id_filters:
                            if pattern.match(resource_id) is not None:
                                id_pass = True
                    else:
                        id_pass = True
                    if resource_type_filters:
                        for resource_type in resource_type_filters:
                            if (
                                EC2_PREFIX_TO_RESOURCE[get_prefix(resource_id)]
                                == resource_type
                            ):
                                type_pass = True
                    else:
                        type_pass = True
                    if value_filters:
                        for pattern in value_filters:
                            if pattern.match(value) is not None:
                                value_pass = True
                    else:
                        value_pass = True
                    if key_pass and id_pass and type_pass and value_pass:
                        add_result = True
                        # If we're not filtering, or we are filtering and this
                if add_result:
                    result = {
                        "resource_id": resource_id,
                        "key": key,
                        "value": value,
                        "resource_type": EC2_PREFIX_TO_RESOURCE.get(
                            get_prefix(resource_id), ""
                        ),
                    }
                    results.append(result)
        return results
