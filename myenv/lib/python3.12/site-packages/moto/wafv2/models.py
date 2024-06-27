import re
from collections import OrderedDict
from typing import TYPE_CHECKING, Any, Dict, List, Optional

from moto.core.base_backend import BackendDict, BaseBackend
from moto.core.common_models import BaseModel
from moto.core.utils import iso_8601_datetime_with_milliseconds
from moto.moto_api._internal import mock_random
from moto.utilities.tagging_service import TaggingService
from moto.utilities.utils import ARN_PARTITION_REGEX, PARTITION_NAMES

from .exceptions import (
    WAFNonexistentItemException,
    WAFOptimisticLockException,
    WAFV2DuplicateItemException,
)
from .utils import make_arn_for_ip_set, make_arn_for_wacl

if TYPE_CHECKING:
    from moto.apigateway.models import Stage


US_EAST_1_REGION = "us-east-1"
GLOBAL_REGION = "global"
APIGATEWAY_REGEX = (
    ARN_PARTITION_REGEX
    + r":apigateway:[a-zA-Z0-9-]+::/restapis/[a-zA-Z0-9]+/stages/[a-zA-Z0-9]+"
)


# TODO: Add remaining properties
class FakeWebACL(BaseModel):
    """
    https://docs.aws.amazon.com/waf/latest/APIReference/API_WebACL.html
    """

    def __init__(
        self,
        name: str,
        arn: str,
        wacl_id: str,
        visibility_config: Dict[str, Any],
        default_action: Dict[str, Any],
        description: Optional[str],
        rules: List[Dict[str, Any]],
    ):
        self.name = name
        self.created_time = iso_8601_datetime_with_milliseconds()
        self.id = wacl_id
        self.arn = arn
        self.description = description or ""
        self.capacity = 3
        self.rules = rules
        self.visibility_config = visibility_config
        self.default_action = default_action
        self.lock_token = str(mock_random.uuid4())[0:6]

    def update(
        self,
        default_action: Optional[Dict[str, Any]],
        rules: Optional[List[Dict[str, Any]]],
        description: Optional[str],
        visibility_config: Optional[Dict[str, Any]],
    ) -> None:
        if default_action is not None:
            self.default_action = default_action
        if rules is not None:
            self.rules = rules
        if description is not None:
            self.description = description
        if visibility_config is not None:
            self.visibility_config = visibility_config
        self.lock_token = str(mock_random.uuid4())[0:6]

    def to_dict(self) -> Dict[str, Any]:
        # Format for summary https://docs.aws.amazon.com/waf/latest/APIReference/API_CreateWebACL.html (response syntax section)
        return {
            "ARN": self.arn,
            "Description": self.description,
            "Id": self.id,
            "Name": self.name,
            "Rules": self.rules,
            "DefaultAction": self.default_action,
            "VisibilityConfig": self.visibility_config,
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
        addresses: List[str],
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

    def update(self, description: Optional[str], addresses: List[str]) -> None:
        if description is not None:
            self.description = description
        self.addresses = addresses

        self.lock_token = str(mock_random.uuid4())[0:6]

    def to_dict(self) -> Dict[str, Any]:
        return {
            "Name": self.name,
            "Id": self.ip_set_id,
            "ARN": self.arn,
            "Description": self.description,
            "IPAddressVersion": self.ip_address_version,
            "Addresses": self.addresses,
            "LockToken": self.lock_token,
        }


class FakeLoggingConfiguration:
    def __init__(
        self,
        arn: str,
        log_destination_configs: List[str],
        redacted_fields: Optional[Dict[str, Any]],
        managed_gy_firewall_manager: Optional[bool],
        logging_filter: Optional[Dict[str, Any]],
    ):
        self.arn = arn
        self.log_destination_configs = log_destination_configs
        self.redacted_fields = redacted_fields
        self.managed_by_firewall_manager = managed_gy_firewall_manager or False
        self.logging_filter = logging_filter

    def to_dict(self) -> Dict[str, Any]:
        return {
            "ResourceArn": self.arn,
            "LogDestinationConfigs": self.log_destination_configs,
            "RedactedFields": self.redacted_fields,
            "ManagedByFirewallManager": self.managed_by_firewall_manager,
            "LoggingFilter": self.logging_filter,
        }


class WAFV2Backend(BaseBackend):
    """
    https://docs.aws.amazon.com/waf/latest/APIReference/API_Operations_AWS_WAFV2.html
    """

    def __init__(self, region_name: str, account_id: str):
        super().__init__(region_name, account_id)
        self.wacls: Dict[str, FakeWebACL] = OrderedDict()
        self.ip_sets: Dict[str, FakeIPSet] = OrderedDict()
        self.logging_configurations: Dict[str, FakeLoggingConfiguration] = OrderedDict()
        self.tagging_service = TaggingService()
        # TODO: self.load_balancers = OrderedDict()

    def associate_web_acl(self, web_acl_arn: str, resource_arn: str) -> None:
        """
        Only APIGateway Stages can be associated at the moment.
        """
        if web_acl_arn not in self.wacls:
            raise WAFNonexistentItemException
        stage = self._find_apigw_stage(resource_arn)
        if stage:
            stage.web_acl_arn = web_acl_arn

    def disassociate_web_acl(self, resource_arn: str) -> None:
        stage = self._find_apigw_stage(resource_arn)
        if stage:
            stage.web_acl_arn = None

    def get_web_acl_for_resource(self, resource_arn: str) -> Optional[FakeWebACL]:
        stage = self._find_apigw_stage(resource_arn)
        if stage and stage.web_acl_arn is not None:
            wacl_arn = stage.web_acl_arn
            return self.wacls.get(wacl_arn)
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

    def create_web_acl(
        self,
        name: str,
        visibility_config: Dict[str, Any],
        default_action: Dict[str, Any],
        scope: str,
        description: str,
        tags: List[Dict[str, str]],
        rules: List[Dict[str, Any]],
    ) -> FakeWebACL:
        """
        The following parameters are not yet implemented: CustomResponseBodies, CaptchaConfig
        """
        wacl_id = str(mock_random.uuid4())
        arn = make_arn_for_wacl(
            name=name,
            account_id=self.account_id,
            region_name=self.region_name,
            wacl_id=wacl_id,
            scope=scope,
        )
        if arn in self.wacls or self._is_duplicate_name(name):
            raise WAFV2DuplicateItemException()
        new_wacl = FakeWebACL(
            name, arn, wacl_id, visibility_config, default_action, description, rules
        )
        self.wacls[arn] = new_wacl
        self.tag_resource(arn, tags)
        return new_wacl

    def delete_web_acl(self, name: str, _id: str) -> None:
        """
        The LockToken-parameter is not yet implemented
        """
        self.wacls = {
            arn: wacl
            for arn, wacl in self.wacls.items()
            if wacl.name != name and wacl.id != _id
        }

    def get_web_acl(self, name: str, _id: str) -> FakeWebACL:
        for wacl in self.wacls.values():
            if wacl.name == name and wacl.id == _id:
                return wacl
        raise WAFNonexistentItemException

    def list_web_acls(self) -> List[Dict[str, Any]]:
        return [wacl.to_dict() for wacl in self.wacls.values()]

    def _is_duplicate_name(self, name: str) -> bool:
        allWaclNames = set(wacl.name for wacl in self.wacls.values())
        return name in allWaclNames

    def list_rule_groups(self) -> List[Any]:
        return []

    def list_tags_for_resource(self, arn: str) -> List[Dict[str, str]]:
        """
        Pagination is not yet implemented
        """
        return self.tagging_service.list_tags_for_resource(arn)["Tags"]

    def tag_resource(self, arn: str, tags: List[Dict[str, str]]) -> None:
        self.tagging_service.tag_resource(arn, tags)

    def untag_resource(self, arn: str, tag_keys: List[str]) -> None:
        self.tagging_service.untag_resource_using_names(arn, tag_keys)

    def update_web_acl(
        self,
        name: str,
        _id: str,
        default_action: Optional[Dict[str, Any]],
        rules: Optional[List[Dict[str, Any]]],
        description: Optional[str],
        visibility_config: Optional[Dict[str, Any]],
    ) -> str:
        """
        The following parameters are not yet implemented: LockToken, CustomResponseBodies, CaptchaConfig
        """
        acl = self.get_web_acl(name, _id)
        acl.update(default_action, rules, description, visibility_config)
        return acl.lock_token

    def create_ip_set(
        self,
        name: str,
        scope: str,
        description: str,
        ip_address_version: str,
        addresses: List[str],
        tags: List[Dict[str, str]],
    ) -> FakeIPSet:
        ip_set_id = str(mock_random.uuid4())
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

    def list_ip_sets(self, scope: str) -> List[FakeIPSet]:
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
        addresses: List[str],
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
        log_destination_configs: List[str],
        redacted_fields: Optional[Dict[str, Any]],
        managed_gy_firewall_manager: bool,
        logging_filter: Dict[str, Any],
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

    def list_logging_configurations(self, scope: str) -> List[FakeLoggingConfiguration]:
        if scope == "CLOUDFRONT":
            scope = "global"
        else:
            scope = self.region_name

        return [
            logging_configuration
            for arn, logging_configuration in self.logging_configurations.items()
            if f":{scope}:" in arn
        ]


wafv2_backends = BackendDict(
    WAFV2Backend, "waf-regional", additional_regions=PARTITION_NAMES
)
