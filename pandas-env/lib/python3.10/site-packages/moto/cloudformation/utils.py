import os
import string
from typing import Any, List
from urllib.parse import urlparse

import yaml

from moto.moto_api._internal import mock_random as random
from moto.utilities.utils import get_partition

from .exceptions import ValidationError


def generate_stack_id(stack_name: str, region: str, account: str) -> str:
    random_id = random.uuid4()
    return f"arn:{get_partition(region)}:cloudformation:{region}:{account}:stack/{stack_name}/{random_id}"


def generate_changeset_id(
    changeset_name: str, region_name: str, account_id: str
) -> str:
    random_id = random.uuid4()
    return f"arn:{get_partition(region_name)}:cloudformation:{region_name}:{account_id}:changeSet/{changeset_name}/{random_id}"


def generate_stackset_id(stackset_name: str) -> str:
    random_id = random.uuid4()
    return f"{stackset_name}:{random_id}"


def generate_stackset_arn(stackset_id: str, region_name: str, account_id: str) -> str:
    return f"arn:{get_partition(region_name)}:cloudformation:{region_name}:{account_id}:stackset/{stackset_id}"


def random_suffix() -> str:
    size = 12
    chars = list(range(10)) + list(string.ascii_uppercase)
    return "".join(str(random.choice(chars)) for x in range(size))


def yaml_tag_constructor(loader: Any, tag: Any, node: Any) -> Any:
    """convert shorthand intrinsic function to full name"""

    def _f(loader: Any, tag: Any, node: Any) -> Any:
        if tag == "!GetAtt":
            if isinstance(node.value, list):
                return node.value
            return node.value.split(".")
        elif type(node) == yaml.SequenceNode:
            return loader.construct_sequence(node)
        else:
            return node.value

    if tag == "!Ref":
        key = "Ref"
    else:
        key = f"Fn::{tag[1:]}"

    return {key: _f(loader, tag, node)}


def validate_template_cfn_lint(template: str) -> List[Any]:
    # Importing cfnlint adds a significant overhead, so we keep it local

    try:
        # Compatibility for cfn-lint 0.x
        # Fail fast with `cfnlint.core.configure_logging` which is removed in cfn-lint 1.x

        from cfnlint.core import configure_logging, get_rules, run_checks
        from cfnlint.decode import decode

        # Save the template to a temporary file -- cfn-lint requires a file
        filename = "file.tmp"
        with open(filename, "w") as file:
            file.write(template)
        abs_filename = os.path.abspath(filename)

        # decode handles both yaml and json
        try:
            template, matches = decode(abs_filename, False)
        except TypeError:
            # As of cfn-lint 0.39.0, the second argument (ignore_bad_template) was dropped
            # https://github.com/aws-cloudformation/cfn-python-lint/pull/1580
            template, matches = decode(abs_filename)

        # Set cfn-lint to info
        configure_logging(None)

        # Initialize the ruleset to be applied (no overrules, no excludes)
        rules = get_rules([], [], [])

        # Use us-east-1 region (spec file) for validation
        regions = ["us-east-1"]

        # Process all the rules and gather the errors
        return run_checks(abs_filename, template, rules, regions)

    except ImportError:
        # Compatibility for cfn-lint 1.x

        from cfnlint.api import lint
        from cfnlint.config import configure_logging
        from cfnlint.core import get_rules

        # Set cfn-lint to info
        configure_logging(None, False)

        # Initialize the ruleset to be applied (no overrules, no excludes)
        rules = get_rules([], [], [])

        # Use us-east-1 region (spec file) for validation
        regions = ["us-east-1"]

        # Process all the rules and gather the errors
        return lint(template, rules, regions)


def get_stack_from_s3_url(template_url: str, account_id: str, partition: str) -> str:
    from moto.s3.models import s3_backends

    template_url_parts = urlparse(template_url)
    if "localhost" in template_url:
        bucket_name, key_name = template_url_parts.path.lstrip("/").split("/", 1)
    else:
        if template_url_parts.netloc.endswith(
            "amazonaws.com"
        ) and template_url_parts.netloc.startswith("s3"):
            # Handle when S3 url uses amazon url with bucket in path
            # Also handles getting region as technically s3 is region'd

            # region = template_url.netloc.split('.')[1]
            bucket_name, key_name = template_url_parts.path.lstrip("/").split("/", 1)
        else:
            bucket_name = template_url_parts.netloc.split(".")[0]
            key_name = template_url_parts.path.lstrip("/")

    key = s3_backends[account_id][partition].get_object(bucket_name, key_name)
    return key.value.decode("utf-8")  # type: ignore[union-attr]


def validate_create_change_set(change_set_name: str) -> None:
    if not (change_set_name and change_set_name[0].isalpha()):
        raise ValidationError(f"Invalid change set name: {change_set_name}")

    if not all(c.isalnum() or c == "-" for c in change_set_name):
        raise ValidationError(f"Invalid change set name: {change_set_name}")

    if len(change_set_name) > 128:
        raise ValidationError(
            f"Change set name exceeds 128 characters: {change_set_name}"
        )

    # Additional validations can be added here later
    return
