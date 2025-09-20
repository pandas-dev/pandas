"""Tag functionality contained in class TaggingService."""

import re
from typing import Dict, List, Optional


class TaggingService:
    """Functionality related to tags, i.e., adding, deleting, testing."""

    def __init__(
        self, tag_name: str = "Tags", key_name: str = "Key", value_name: str = "Value"
    ):
        self.tag_name = tag_name
        self.key_name = key_name
        self.value_name = value_name
        self.tags: Dict[str, Dict[str, Optional[str]]] = {}

    def get_tag_dict_for_resource(self, arn: str) -> Dict[str, str]:
        """Return dict of key/value pairs vs. list of key/values dicts."""
        result = {}
        if self.has_tags(arn):
            for key, val in self.tags[arn].items():
                result[key] = val
        return result  # type: ignore

    def list_tags_for_resource(self, arn: str) -> Dict[str, List[Dict[str, str]]]:
        """Return list of tags inside dict with key of "tag_name".

        Useful for describe functions; this return value can be added to
        dictionary returned from a describe function.
        """
        result = []
        if self.has_tags(arn):
            for key, val in self.tags[arn].items():
                result.append({self.key_name: key, self.value_name: val})
        return {self.tag_name: result}  # type: ignore

    def delete_all_tags_for_resource(self, arn: str) -> None:
        """Delete all tags associated with given ARN."""
        if self.has_tags(arn):
            del self.tags[arn]

    def has_tags(self, arn: str) -> bool:
        """Return True if the ARN has any associated tags, False otherwise."""
        return arn in self.tags

    def tag_resource(self, arn: str, tags: Optional[List[Dict[str, str]]]) -> None:
        """Store associated list of dicts with ARN.

        Note: the storage is internal to this class instance.
        """
        if not tags:
            return
        if arn not in self.tags:
            self.tags[arn] = {}
        for tag in tags:
            if self.value_name in tag:
                self.tags[arn][tag[self.key_name]] = tag[self.value_name]
            else:
                self.tags[arn][tag[self.key_name]] = None

    def copy_tags(self, from_arn: str, to_arn: str) -> None:
        """Copy stored list of tags associated with one ARN to another ARN.

        Note: the storage is internal to this class instance.
        """
        if self.has_tags(from_arn):
            self.tag_resource(
                to_arn, self.list_tags_for_resource(from_arn)[self.tag_name]
            )

    def untag_resource_using_names(self, arn: str, tag_names: List[str]) -> None:
        """Remove tags associated with ARN using key names in 'tag_names'."""
        for name in tag_names:
            if name in self.tags.get(arn, {}):
                del self.tags[arn][name]

    def untag_resource_using_tags(self, arn: str, tags: List[Dict[str, str]]) -> None:
        """Remove tags associated with ARN using key/value pairs in 'tags'."""
        current_tags = self.tags.get(arn, {})
        for tag in tags:
            if self.key_name in tag:
                if tag[self.key_name] in current_tags:
                    if self.value_name in tag:
                        if current_tags[tag[self.key_name]] != tag[self.value_name]:
                            continue
                    # If both key and value are provided, match both before deletion
                    del current_tags[tag[self.key_name]]

    def extract_tag_names(self, tags: List[Dict[str, str]]) -> List[str]:
        """Return list of key names in list of 'tags' key/value dicts."""
        results: List[str] = []
        if len(tags) == 0:
            return results
        for tag in tags:
            if self.key_name in tag:
                results.append(tag[self.key_name])
        return results

    def flatten_tag_list(self, tags: List[Dict[str, str]]) -> Dict[str, Optional[str]]:
        """Return dict of key/value pairs with 'tag_name', 'value_name'."""
        result: Dict[str, Optional[str]] = {}
        for tag in tags:
            if self.value_name in tag:
                result[tag[self.key_name]] = tag[self.value_name]
            else:
                result[tag[self.key_name]] = None
        return result

    def validate_tags(self, tags: List[Dict[str, str]], limit: int = 0) -> str:
        """Returns error message if tags in 'tags' list of dicts are invalid.

        The validation does not include a check for duplicate keys.
        Duplicate keys are not always an error and the error message isn't
        consistent across services, so this should be a separate check.

        If limit is provided, then the number of tags will be checked.
        """
        errors = []
        key_regex = re.compile(r"^(?!aws:)([\w\s\d_.:/=+\-@]*)$")
        value_regex = re.compile(r"^([\w\s\d_.:/=+\-@]*)$")

        # AWS only outputs one error for all keys and one for all values.
        for idx, tag in enumerate(tags, 1):
            for tag_key, tag_value in tag.items():
                if tag_key == self.key_name:
                    # Validation for len(tag_key) >= 1 is done by botocore.
                    if len(tag_value) > 128:
                        errors.append(
                            f"Value '{tag_value}' at 'tags.{idx}.member.key' "
                            f"failed to satisfy constraint: Member must have "
                            f"length less than or equal to 128"
                        )
                    if not re.match(key_regex, tag_value):
                        errors.append(
                            f"Value '{tag_value}' at 'tags.{idx}.member.key' "
                            f"failed to satisfy constraint: Member must "
                            f"satisfy regular expression pattern: "
                            r"^(?!aws:)[{a-zA-Z0-9 }_.://=+-@%]*$"
                        )
                elif tag_key == self.value_name:
                    # Validation for len(tag_value) >= 0 is nonsensical.
                    if len(tag_value) > 256:
                        errors.append(
                            f"Value '{tag_value}' at 'tags.{idx}.member.value' "
                            f"failed to satisfy constraint: Member must have "
                            f"length less than or equal to 256"
                            # Member must have length greater than or equal to 0, "
                        )
                    if not re.match(value_regex, tag_value):
                        errors.append(
                            f"Value '{tag_value}' at 'tags.{idx}.member.value' "
                            f"failed to satisfy constraint: Member must satisfy "
                            f"regular expression pattern: "
                            r"^[{a-zA-Z0-9 }_.://=+-@%]*$"
                        )

        if limit and len(tags) > limit:
            errors.append(
                f"Value '{tags}' at 'tags' failed to satisfy constraint: "
                f"Member must have length less than or equal to {limit}"
            )

        errors_len = len(errors)
        return (
            (
                f"{errors_len} validation error{'s' if len(errors) > 1 else ''} "
                f"detected: {'; '.join(errors)}"
            )
            if errors
            else ""
        )

    @staticmethod
    def convert_dict_to_tags_input(
        tags: Optional[Dict[str, str]],
    ) -> List[Dict[str, str]]:
        """Given a dictionary, return generic boto params for tags"""
        if not tags:
            return []
        return [{"Key": k, "Value": v} for (k, v) in tags.items()]
