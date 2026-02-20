import re
from collections import OrderedDict
from typing import TYPE_CHECKING, Any, Optional, Union

from moto.core.base_backend import BackendDict, BaseBackend
from moto.core.common_models import BaseModel
from moto.core.utils import (
    camelcase_to_underscores,
    iso_8601_datetime_with_milliseconds,
)
from moto.moto_api._internal import mock_random
from moto.utilities.paginator import paginate
from moto.utilities.tagging_service import TaggingService
from moto.utilities.utils import ARN_PARTITION_REGEX, PARTITION_NAMES

from .exceptions import (
    WAFNonexistentItemException,
    WAFOptimisticLockException,
    WAFV2DuplicateItemException,
    WAFV2InsufficientInformationException,
)
from .utils import (
    make_arn_for_ip_set,
    make_arn_for_regex_pattern_set,
    make_arn_for_rule_group,
    make_arn_for_wacl,
)

if TYPE_CHECKING:
    from moto.apigateway.models import Stage


US_EAST_1_REGION = "us-east-1"
GLOBAL_REGION = "global"
APIGATEWAY_REGEX = (
    ARN_PARTITION_REGEX
    + r":apigateway:[a-zA-Z0-9-]+::/restapis/[a-zA-Z0-9]+/stages/[a-zA-Z0-9]+"
)

PAGINATION_MODEL = {
    "list_ip_sets": {
        "input_token": "next_marker",
        "limit_key": "limit",
        "limit_default": 100,
        "unique_attribute": "arn",
    },
    "list_logging_configurations": {
        "input_token": "next_marker",
        "limit_key": "limit",
        "limit_default": 100,
        "unique_attribute": "arn",
    },
    "list_regex_pattern_sets": {
        "input_token": "next_marker",
        "limit_key": "limit",
        "limit_default": 100,
        "unique_attribute": "arn",
    },
    "list_rule_groups": {
        "input_token": "next_marker",
        "limit_key": "limit",
        "limit_default": 100,
        "unique_attribute": "arn",
    },
    "list_tags_for_resource": {
        "input_token": "next_marker",
        "limit_key": "limit",
        "limit_default": 100,
        "unique_attribute": "Key",
    },
    "list_web_acls": {
        "input_token": "next_marker",
        "limit_key": "limit",
        "limit_default": 100,
        "unique_attribute": "arn",
    },
}


class FakeRule(BaseModel):
    def __init__(
        self,
        name: str,
        priority: int,
        statement: dict[str, Any],
        visibility_config: dict[str, Union[str, bool]],
        action: Optional[dict[str, Any]] = None,
        captcha_config: Optional[dict[str, dict[str, int]]] = None,
        challenge_config: Optional[dict[str, dict[str, int]]] = None,
        override_action: Optional[dict[str, Any]] = None,
        rule_labels: Optional[list[dict[str, str]]] = None,
    ):
        self.name = name
        self.priority = priority
        self.statement = statement
        self.visibility_config = visibility_config
        self.action = action
        self.captcha_config = captcha_config
        self.challenge_config = challenge_config
        self.override_action = override_action
        self.rule_labels = rule_labels

    def get_consumed_label(self) -> Optional[str]:
        return self.statement.get("LabelMatchStatement", {}).get("Key")

    def get_available_labels(self) -> Optional[list[str]]:
        return [r["Name"] for r in self.rule_labels] if self.rule_labels else None

    def to_dict(self) -> dict[str, Any]:
        return {
            "Name": self.name,
            "Action": self.action,
            "Priority": self.priority,
            "Statement": self.statement,
            "VisibilityConfig": self.visibility_config,
            "CaptchaConfig": self.captcha_config,
            "ChallengeConfig": self.challenge_config,
            "OverrideAction": self.override_action,
            "RuleLabels": self.rule_labels,
        }


class FakeWebACL(BaseModel):
    """
    https://docs.aws.amazon.com/waf/latest/APIReference/API_WebACL.html
    """

    def __init__(
        self,
        name: str,
        account: str,
        arn: str,
        wacl_id: str,
        visibility_config: dict[str, Any],
        default_action: dict[str, Any],
        description: Optional[str],
        rules: Optional[list[FakeRule]],
        association_config: Optional[dict[str, Any]] = None,
        captcha_config: Optional[dict[str, Any]] = None,
        challenge_config: Optional[dict[str, Any]] = None,
        custom_response_bodies: Optional[dict[str, Any]] = None,
        token_domains: Optional[list[str]] = None,
    ):
        self.name = name
        self.account = account
        self.created_time = iso_8601_datetime_with_milliseconds()
        self.id = wacl_id
        self.arn = arn
        self.description = description or ""
        self.capacity = 3
        self.rules = rules or []
        self.visibility_config = visibility_config
        self.default_action = default_action
        self.lock_token = self._generate_lock_token()
        self.associated_resources: list[str] = []
        self.association_config = association_config
        self.captcha_config = captcha_config
        self.challenge_config = challenge_config
        self.custom_response_bodies = custom_response_bodies
        self.token_domains = token_domains
        self.label_namespace = self._get_label_namespace()

    def _generate_lock_token(self) -> str:
        return str(mock_random.uuid4())

    def _get_label_namespace(self) -> str:
        return f"awswaf:{self.account}:webacl:{self.name}:"

    def update(
        self,
        default_action: Optional[dict[str, Any]],
        rules: Optional[list[FakeRule]],
        description: Optional[str],
        visibility_config: Optional[dict[str, Any]],
        custom_response_bodies: Optional[dict[str, Any]],
        captcha_config: Optional[dict[str, Any]],
        challenge_config: Optional[dict[str, Any]],
        token_domains: Optional[list[str]],
        association_config: Optional[dict[str, Any]],
    ) -> None:
        if default_action is not None:
            self.default_action = default_action
        if rules is not None:
            self.rules = rules
        if description is not None:
            self.description = description
        if visibility_config is not None:
            self.visibility_config = visibility_config
        if custom_response_bodies is not None:
            self.custom_response_bodies = custom_response_bodies
        if captcha_config is not None:
            self.captcha_config = captcha_config
        if challenge_config is not None:
            self.challenge_config = challenge_config
        if token_domains is not None:
            self.token_domains = token_domains
        if association_config is not None:
            self.association_config = association_config
        self.lock_token = self._generate_lock_token()

    def to_short_dict(self) -> dict[str, Any]:
        # Format for summary https://docs.aws.amazon.com/waf/latest/APIReference/API_CreateWebACL.html (response syntax section)
        return {
            "ARN": self.arn,
            "Description": self.description,
            "Id": self.id,
            "Name": self.name,
            "LockToken": self.lock_token,
        }

    def to_dict(self) -> dict[str, Any]:
        return {
            "Name": self.name,
            "Id": self.id,
            "ARN": self.arn,
            "Description": self.description,
            "Rules": [r.to_dict() for r in self.rules],
            "DefaultAction": self.default_action,
            "VisibilityConfig": self.visibility_config,
            "Capacity": self.capacity,
            "LockToken": self.lock_token,
            "AssociationConfig": self.association_config,
            "CaptchaConfig": self.captcha_config,
            "ChallengeConfig": self.challenge_config,
            "CustomResponseBodies": self.custom_response_bodies,
            "TokenDomains": self.token_domains,
            "LabelNamespace": self.label_namespace,
        }


class FakeIPSet(BaseModel):
    """
    https://docs.aws.amazon.com/waf/latest/APIReference/API_IPSet.html
    """

    def __init__(
        self,
        arn: str,
        ip_set_id: str,
        ip_address_version: str,
        addresses: list[str],
        name: str,
        description: str,
        scope: str,
    ):
        self.name = name
        self.ip_set_id = ip_set_id
        self.arn = arn
        self.addresses = addresses
        self.description = description
        self.ip_address_version = ip_address_version
        self.scope = scope

        self.lock_token = str(mock_random.uuid4())[0:6]

    def update(self, description: Optional[str], addresses: list[str]) -> None:
        if description is not None:
            self.description = description
        self.addresses = addresses

        self.lock_token = str(mock_random.uuid4())[0:6]

    def to_dict(self) -> dict[str, Any]:
        return {
            "Name": self.name,
            "Id": self.ip_set_id,
            "ARN": self.arn,
            "Description": self.description,
            "IPAddressVersion": self.ip_address_version,
            "Addresses": self.addresses,
            "LockToken": self.lock_token,
        }


class FakeLoggingConfiguration(BaseModel):
    def __init__(
        self,
        arn: str,
        log_destination_configs: list[str],
        redacted_fields: Optional[dict[str, Any]],
        managed_gy_firewall_manager: Optional[bool],
        logging_filter: Optional[dict[str, Any]],
    ):
        self.arn = arn
        self.log_destination_configs = log_destination_configs
        self.redacted_fields = redacted_fields
        self.managed_by_firewall_manager = managed_gy_firewall_manager or False
        self.logging_filter = logging_filter

    def to_dict(self) -> dict[str, Any]:
        return {
            "ResourceArn": self.arn,
            "LogDestinationConfigs": self.log_destination_configs,
            "RedactedFields": self.redacted_fields,
            "ManagedByFirewallManager": self.managed_by_firewall_manager,
            "LoggingFilter": self.logging_filter,
        }


class FakeRuleGroup(BaseModel):
    def __init__(
        self,
        account: str,
        name: str,
        id: str,
        arn: str,
        scope: str,
        capacity: int,
        visibility_config: dict[str, Any],
        description: Optional[str],
        rules: Optional[list[FakeRule]],
        custom_response_bodies: Optional[dict[str, Any]] = None,
    ):
        self.account = account
        self.name = name
        self.id = id
        self.arn = arn
        self.lock_token = self._generate_lock_token()
        self.scope = scope
        self.capacity = capacity
        self.description = description or ""
        self.rules = rules or []
        self.visibility_config = visibility_config
        self.custom_response_bodies = custom_response_bodies
        self.label_namespace = self._get_label_namespace()
        self.available_labels = self._get_available_labels()
        self.consumed_labels = self._get_consumed_labels()

    def _generate_lock_token(self) -> str:
        return str(mock_random.uuid4())

    def _get_available_labels(self) -> Optional[list[dict[str, str]]]:
        return (
            [
                {"Name": f"{self.label_namespace}{label}"}
                for rule in self.rules
                for label in (rule.get_available_labels() or [])
            ]
            if self.rules
            else None
        )

    def _get_consumed_labels(self) -> Optional[list[dict[str, str]]]:
        return (
            [
                {"Name": f"{self.label_namespace}{label}"}
                for rule in self.rules
                if (label := rule.get_consumed_label())
            ]
            if self.rules
            else None
        )

    def _get_label_namespace(self) -> str:
        return f"awswaf:{self.account}:rulegroup:{self.name}:"

    def update(
        self,
        description: Optional[str],
        rules: Optional[list[FakeRule]],
        visibility_config: Optional[dict[str, Any]],
        custom_response_bodies: Optional[dict[str, Any]],
    ) -> str:
        if description is not None:
            self.description = description
        if rules is not None:
            self.rules = rules
        if visibility_config is not None:
            self.visibility_config = visibility_config
        if custom_response_bodies is not None:
            self.custom_response_bodies = custom_response_bodies

        self.lock_token = self._generate_lock_token()
        return self.lock_token

    def to_short_dict(self) -> dict[str, Any]:
        return {
            "Name": self.name,
            "Id": self.id,
            "Description": self.description,
            "LockToken": self.lock_token,
            "ARN": self.arn,
        }

    def to_dict(self) -> dict[str, Any]:
        return {
            "Name": self.name,
            "Id": self.id,
            "Capacity": self.capacity,
            "ARN": self.arn,
            "Description": self.description,
            "Rules": [r.to_dict() for r in self.rules],
            "VisibilityConfig": self.visibility_config,
            "LabelNamespace": self.label_namespace,
            "AvailableLabels": self.available_labels,
            "ConsumedLabels": self.consumed_labels,
            "CustomResponseBodies": self.custom_response_bodies,
            "LockToken": self.lock_token,
        }


class FakeRegexPatternSet(BaseModel):
    """
    https://docs.aws.amazon.com/waf/latest/APIReference/API_RegexPatternSet.html
    """

    def __init__(
        self,
        name: str,
        account: str,
        id: str,
        arn: str,
        description: str,
        regular_expressions: list[dict[str, str]],
        scope: str,
    ):
        self.name = name
        self.account = account
        self.id = id
        self.arn = arn
        self.description = description
        self.regular_expressions = regular_expressions
        self.scope = scope
        self.lock_token = str(mock_random.uuid4())

    def update(
        self,
        description: Optional[str],
        regular_expressions: Optional[list[dict[str, str]]],
    ) -> None:
        if description is not None:
            self.description = description
        if regular_expressions is not None:
            self.regular_expressions = regular_expressions
        self.lock_token = str(mock_random.uuid4())

    def to_dict(self) -> dict[str, Any]:
        return {
            "Name": self.name,
            "Id": self.id,
            "ARN": self.arn,
            "Description": self.description,
            "RegularExpressionList": self.regular_expressions,
            "LockToken": self.lock_token,
        }

    def to_short_dict(self) -> dict[str, Any]:
        return {
            "Name": self.name,
            "Id": self.id,
            "Description": self.description,
            "ARN": self.arn,
            "LockToken": self.lock_token,
        }


class WAFV2Backend(BaseBackend):
    """
    https://docs.aws.amazon.com/waf/latest/APIReference/API_Operations_AWS_WAFV2.html
    """

    def __init__(self, region_name: str, account_id: str):
        super().__init__(region_name, account_id)
        self.wacls: dict[str, FakeWebACL] = OrderedDict()
        self.ip_sets: dict[str, FakeIPSet] = OrderedDict()
        self.logging_configurations: dict[str, FakeLoggingConfiguration] = OrderedDict()
        self.rule_groups: dict[str, FakeRuleGroup] = OrderedDict()
        self.regex_pattern_sets: dict[str, FakeRegexPatternSet] = OrderedDict()
        self.tagging_service = TaggingService()
        # TODO: self.load_balancers = OrderedDict()

    def associate_web_acl(self, web_acl_arn: str, resource_arn: str) -> None:
        web_acl = self.wacls.get(web_acl_arn)
        if not web_acl:
            raise WAFNonexistentItemException
        web_acl.associated_resources.append(resource_arn)
        # Special Case - APIGateway wants to know about the WebACL it's associated to
        stage = self._find_apigw_stage(resource_arn)
        if stage:
            stage.web_acl_arn = web_acl_arn

    def disassociate_web_acl(self, resource_arn: str) -> None:
        for web_acl in self.wacls.values():
            if resource_arn in web_acl.associated_resources:
                web_acl.associated_resources.remove(resource_arn)
                break
        stage = self._find_apigw_stage(resource_arn)
        if stage:
            stage.web_acl_arn = None

    def get_web_acl_for_resource(self, resource_arn: str) -> Optional[FakeWebACL]:
        for wacl in self.wacls.values():
            if resource_arn in wacl.associated_resources:
                return wacl
        return None

    def _find_apigw_stage(self, resource_arn: str) -> Optional["Stage"]:  # type: ignore
        try:
            if re.search(APIGATEWAY_REGEX, resource_arn):
                region = resource_arn.split(":")[3]
                rest_api_id = resource_arn.split("/")[-3]
                stage_name = resource_arn.split("/")[-1]

                from moto.apigateway import apigateway_backends

                apigw = apigateway_backends[self.account_id][region]
                return apigw.get_stage(rest_api_id, stage_name)
        except:  # noqa: E722 Do not use bare except
            return None

    def _generate_id(self) -> str:
        return str(mock_random.uuid4())

    def create_web_acl(
        self,
        name: str,
        visibility_config: dict[str, Any],
        default_action: dict[str, Any],
        scope: str,
        description: str,
        tags: list[dict[str, str]],
        rules: list[dict[str, Any]],
        association_config: Optional[dict[str, Any]],
        captcha_config: Optional[dict[str, Any]],
        challenge_config: Optional[dict[str, Any]],
        custom_response_bodies: Optional[dict[str, Any]],
        token_domains: Optional[list[str]],
    ) -> FakeWebACL:
        wacl_id = self._generate_id()
        arn = make_arn_for_wacl(
            name=name,
            account_id=self.account_id,
            region_name=self.region_name,
            _id=wacl_id,
            scope=scope,
        )
        if arn in self.wacls or self._is_duplicate_name(name):
            raise WAFV2DuplicateItemException()
        rule_objs = [
            FakeRule(**{camelcase_to_underscores(k): v for k, v in rule.items()})
            for rule in rules
        ]
        new_wacl = FakeWebACL(
            name=name,
            account=self.account_id,
            arn=arn,
            wacl_id=wacl_id,
            visibility_config=visibility_config,
            default_action=default_action,
            description=description,
            rules=rule_objs,
            association_config=association_config,
            captcha_config=captcha_config,
            challenge_config=challenge_config,
            custom_response_bodies=custom_response_bodies,
            token_domains=token_domains,
        )
        self.wacls[arn] = new_wacl
        self.tag_resource(arn, tags)
        return new_wacl

    def delete_web_acl(self, name: str, scope: str, _id: str, lock_token: str) -> None:
        arn = make_arn_for_wacl(
            name=name,
            account_id=self.account_id,
            region_name=self.region_name,
            _id=_id,
            scope=scope,
        )
        if arn not in self.wacls:
            raise WAFNonexistentItemException
        wacl = self.wacls[arn]
        if wacl.lock_token != lock_token:
            raise WAFOptimisticLockException()
        self.wacls.pop(arn)

    def get_web_acl(self, name: str, _id: str) -> FakeWebACL:
        for wacl in self.wacls.values():
            if wacl.name == name and wacl.id == _id:
                return wacl
        raise WAFNonexistentItemException

    @paginate(PAGINATION_MODEL)  # type: ignore
    def list_web_acls(self) -> list[FakeWebACL]:
        return list(self.wacls.values())

    def _is_duplicate_name(self, name: str) -> bool:
        all_wacl_names = {wacl.name for wacl in self.wacls.values()}
        return name in all_wacl_names

    @paginate(PAGINATION_MODEL)  # type: ignore
    def list_rule_groups(self, scope: str) -> list[FakeRuleGroup]:
        rule_groups = [
            group for group in self.rule_groups.values() if group.scope == scope
        ]
        return rule_groups

    @paginate(PAGINATION_MODEL)  # type: ignore
    def list_tags_for_resource(self, arn: str) -> list[dict[str, str]]:
        return self.tagging_service.list_tags_for_resource(arn)["Tags"]

    def tag_resource(self, arn: str, tags: list[dict[str, str]]) -> None:
        self.tagging_service.tag_resource(arn, tags)

    def untag_resource(self, arn: str, tag_keys: list[str]) -> None:
        self.tagging_service.untag_resource_using_names(arn, tag_keys)

    def update_web_acl(
        self,
        name: str,
        _id: str,
        default_action: Optional[dict[str, Any]],
        rules: Optional[list[dict[str, Any]]],
        description: Optional[str],
        visibility_config: Optional[dict[str, Any]],
        lock_token: str,
        custom_response_bodies: Optional[dict[str, Any]],
        captcha_config: Optional[dict[str, Any]],
        challenge_config: Optional[dict[str, Any]],
        token_domains: Optional[list[str]],
        association_config: Optional[dict[str, Any]],
    ) -> str:
        acl = self.get_web_acl(name, _id)
        if acl.lock_token != lock_token:
            raise WAFOptimisticLockException()
        rule_objs = (
            [
                FakeRule(**{camelcase_to_underscores(k): v for k, v in rule.items()})
                for rule in rules
            ]
            if rules
            else None
        )
        acl.update(
            default_action,
            rule_objs,
            description,
            visibility_config,
            custom_response_bodies,
            captcha_config,
            challenge_config,
            token_domains,
            association_config,
        )
        return acl.lock_token

    def create_ip_set(
        self,
        name: str,
        scope: str,
        description: str,
        ip_address_version: str,
        addresses: list[str],
        tags: list[dict[str, str]],
    ) -> FakeIPSet:
        ip_set_id = self._generate_id()
        arn = make_arn_for_ip_set(
            name=name,
            account_id=self.account_id,
            region_name=self.region_name,
            _id=ip_set_id,
            scope=scope,
        )

        new_ip_set = FakeIPSet(
            arn,
            ip_set_id,
            ip_address_version,
            addresses,
            name,
            description,
            scope,
        )
        self.ip_sets[arn] = new_ip_set
        self.tag_resource(arn, tags)
        return new_ip_set

    def delete_ip_set(self, name: str, scope: str, _id: str, lock_token: str) -> None:
        arn = make_arn_for_ip_set(
            name=name,
            account_id=self.account_id,
            region_name=self.region_name,
            _id=_id,
            scope=scope,
        )

        if arn not in self.ip_sets:
            raise WAFNonexistentItemException()

        if lock_token != self.ip_sets[arn].lock_token:
            raise WAFOptimisticLockException()

        self.ip_sets.pop(arn)

    @paginate(PAGINATION_MODEL)  # type: ignore
    def list_ip_sets(self, scope: str) -> list[FakeIPSet]:
        ip_sets = [
            ip_set for arn, ip_set in self.ip_sets.items() if ip_set.scope == scope
        ]
        return ip_sets

    def get_ip_set(self, name: str, scope: str, _id: str) -> FakeIPSet:
        arn = make_arn_for_ip_set(
            name=name,
            account_id=self.account_id,
            region_name=self.region_name,
            _id=_id,
            scope=scope,
        )
        if arn not in self.ip_sets:
            raise WAFNonexistentItemException()

        return self.ip_sets[arn]

    def update_ip_set(
        self,
        name: str,
        scope: str,
        _id: str,
        description: Optional[str],
        addresses: list[str],
        lock_token: str,
    ) -> FakeIPSet:
        arn = make_arn_for_ip_set(
            name=name,
            account_id=self.account_id,
            region_name=self.region_name,
            _id=_id,
            scope=scope,
        )

        if not (ip_set := self.ip_sets.get(arn)):
            raise WAFNonexistentItemException()

        if ip_set.lock_token != lock_token:
            raise WAFOptimisticLockException()

        ip_set.update(description, addresses)

        return ip_set

    def put_logging_configuration(
        self,
        arn: str,
        log_destination_configs: list[str],
        redacted_fields: Optional[dict[str, Any]],
        managed_gy_firewall_manager: bool,
        logging_filter: dict[str, Any],
    ) -> FakeLoggingConfiguration:
        logging_configuration = FakeLoggingConfiguration(
            arn,
            log_destination_configs,
            redacted_fields,
            managed_gy_firewall_manager,
            logging_filter,
        )
        self.logging_configurations[arn] = logging_configuration
        return logging_configuration

    def delete_logging_configuration(self, arn: str) -> None:
        if not self.logging_configurations.get(arn):
            raise WAFNonexistentItemException()
        self.logging_configurations.pop(arn)

    def get_logging_configuration(self, arn: str) -> FakeLoggingConfiguration:
        if not (logging_configuration := self.logging_configurations.get(arn)):
            raise WAFNonexistentItemException()
        return logging_configuration

    @paginate(PAGINATION_MODEL)  # type: ignore
    def list_logging_configurations(self, scope: str) -> list[FakeLoggingConfiguration]:
        if scope == "CLOUDFRONT":
            scope = "global"
        else:
            scope = self.region_name

        return [
            logging_configuration
            for arn, logging_configuration in self.logging_configurations.items()
            if f":{scope}:" in arn
        ]

    def create_rule_group(
        self,
        name: str,
        scope: str,
        capacity: int,
        description: Optional[str],
        rules: Optional[list[dict[str, Any]]],
        visibility_config: dict[str, Union[bool, str]],
        tags: Optional[list[dict[str, str]]],
        custom_response_bodies: Optional[dict[str, str]],
    ) -> FakeRuleGroup:
        id = self._generate_id()
        arn = make_arn_for_rule_group(
            name, self.account_id, self.region_name, id, scope
        )
        if name in {group.name for group in self.rule_groups.values()}:
            raise WAFV2DuplicateItemException()
        rules_objs = (
            [
                FakeRule(**{camelcase_to_underscores(k): v for k, v in rule.items()})
                for rule in rules
            ]
            if rules
            else None
        )
        rule_group = FakeRuleGroup(
            name=name,
            account=self.account_id,
            id=id,
            arn=arn,
            scope=scope,
            capacity=capacity,
            visibility_config=visibility_config,
            description=description,
            rules=rules_objs,
            custom_response_bodies=custom_response_bodies,
        )
        self.rule_groups[arn] = rule_group
        self.tagging_service.tag_resource(arn, tags)
        return rule_group

    def update_rule_group(
        self,
        name: str,
        scope: str,
        id: str,
        description: Optional[str],
        rules: Optional[list[dict[str, Any]]],
        visibility_config: dict[str, Any],
        lock_token: str,
        custom_response_bodies: Optional[dict[str, Any]],
    ) -> FakeRuleGroup:
        arn = make_arn_for_rule_group(
            name, self.account_id, self.region_name, id, scope
        )
        if not (group := self.rule_groups.get(arn)):
            raise WAFNonexistentItemException()
        if group.lock_token != lock_token:
            raise WAFOptimisticLockException()
        rules_objs = (
            [
                FakeRule(**{camelcase_to_underscores(k): v for k, v in rule.items()})
                for rule in rules
            ]
            if rules
            else None
        )
        group.update(description, rules_objs, visibility_config, custom_response_bodies)

        return group

    def delete_rule_group(
        self, name: str, scope: str, id: str, lock_token: str
    ) -> None:
        arn = make_arn_for_rule_group(
            name, self.account_id, self.region_name, id, scope
        )
        if not (group := self.rule_groups.get(arn)):
            raise WAFNonexistentItemException()
        if group.lock_token != lock_token:
            raise WAFOptimisticLockException()
        self.rule_groups.pop(arn)
        return

    def get_rule_group(
        self,
        name: Optional[str],
        scope: Optional[str],
        id: Optional[str],
        arn: Optional[str],
    ) -> FakeRuleGroup:
        if not arn and not (name and scope and id):
            raise WAFV2InsufficientInformationException(name, scope, id, arn)
        else:
            arn = arn or make_arn_for_rule_group(
                name,  # type: ignore[arg-type]
                self.account_id,
                self.region_name,
                id,  # type: ignore[arg-type]
                scope,  # type: ignore[arg-type]
            )
        if not (group := self.rule_groups.get(arn)):
            raise WAFNonexistentItemException()
        return group

    def _is_duplicate_regex_pattern_set(self, name: str, scope: str) -> bool:
        """Check if a regex pattern set with the same name exists in the same scope"""
        return any(
            pattern_set.name == name and pattern_set.scope == scope
            for pattern_set in self.regex_pattern_sets.values()
        )

    def create_regex_pattern_set(
        self,
        name: str,
        scope: str,
        description: str,
        regular_expressions: list[dict[str, str]],
        tags: list[dict[str, str]],
    ) -> FakeRegexPatternSet:
        regex_pattern_set_id = self._generate_id()
        arn = make_arn_for_regex_pattern_set(
            name=name,
            account_id=self.account_id,
            region_name=self.region_name,
            _id=regex_pattern_set_id,
            scope=scope,
        )

        if self._is_duplicate_regex_pattern_set(name, scope):
            raise WAFV2DuplicateItemException()

        new_regex_pattern_set = FakeRegexPatternSet(
            name=name,
            account=self.account_id,
            id=regex_pattern_set_id,
            arn=arn,
            description=description,
            regular_expressions=regular_expressions,
            scope=scope,
        )
        self.regex_pattern_sets[arn] = new_regex_pattern_set
        self.tag_resource(arn, tags)
        return new_regex_pattern_set

    def get_regex_pattern_set(
        self,
        name: str,
        scope: str,
        id: str,
    ) -> FakeRegexPatternSet:
        arn = make_arn_for_regex_pattern_set(
            name=name,
            account_id=self.account_id,
            region_name=self.region_name,
            _id=id,
            scope=scope,
        )
        if not (pattern_set := self.regex_pattern_sets.get(arn)):
            raise WAFNonexistentItemException()
        return pattern_set

    def update_regex_pattern_set(
        self,
        name: str,
        scope: str,
        id: str,
        description: Optional[str],
        regular_expressions: Optional[list[dict[str, str]]],
        lock_token: str,
    ) -> FakeRegexPatternSet:
        arn = make_arn_for_regex_pattern_set(
            name=name,
            account_id=self.account_id,
            region_name=self.region_name,
            _id=id,
            scope=scope,
        )
        if not (pattern_set := self.regex_pattern_sets.get(arn)):
            raise WAFNonexistentItemException()
        if pattern_set.lock_token != lock_token:
            raise WAFOptimisticLockException()
        pattern_set.update(description, regular_expressions)
        return pattern_set

    def delete_regex_pattern_set(
        self,
        name: str,
        scope: str,
        id: str,
        lock_token: str,
    ) -> None:
        arn = make_arn_for_regex_pattern_set(
            name=name,
            account_id=self.account_id,
            region_name=self.region_name,
            _id=id,
            scope=scope,
        )
        if not (pattern_set := self.regex_pattern_sets.get(arn)):
            raise WAFNonexistentItemException()
        if pattern_set.lock_token != lock_token:
            raise WAFOptimisticLockException()
        self.regex_pattern_sets.pop(arn)

    @paginate(PAGINATION_MODEL)
    def list_regex_pattern_sets(self, scope: str) -> list[FakeRegexPatternSet]:
        return [
            pattern_set
            for pattern_set in self.regex_pattern_sets.values()
            if pattern_set.scope == scope
        ]


wafv2_backends = BackendDict(WAFV2Backend, "wafv2", additional_regions=PARTITION_NAMES)
