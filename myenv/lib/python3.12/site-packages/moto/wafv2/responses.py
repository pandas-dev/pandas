import json

from moto.core.common_types import TYPE_RESPONSE
from moto.core.responses import BaseResponse

from ..moto_api._internal import mock_random
from .models import GLOBAL_REGION, WAFV2Backend, wafv2_backends


class WAFV2Response(BaseResponse):
    def __init__(self) -> None:
        super().__init__(service_name="wafv2")

    @property
    def wafv2_backend(self) -> WAFV2Backend:
        return wafv2_backends[self.current_account][self.region]

    def associate_web_acl(self) -> TYPE_RESPONSE:
        body = json.loads(self.body)
        web_acl_arn = body["WebACLArn"]
        resource_arn = body["ResourceArn"]
        self.wafv2_backend.associate_web_acl(web_acl_arn, resource_arn)
        return 200, {}, "{}"

    def disassociate_web_acl(self) -> TYPE_RESPONSE:
        body = json.loads(self.body)
        resource_arn = body["ResourceArn"]
        self.wafv2_backend.disassociate_web_acl(resource_arn)
        return 200, {}, "{}"

    def get_web_acl_for_resource(self) -> TYPE_RESPONSE:
        body = json.loads(self.body)
        resource_arn = body["ResourceArn"]
        web_acl = self.wafv2_backend.get_web_acl_for_resource(resource_arn)
        response = {"WebACL": web_acl.to_dict() if web_acl else None}
        response_headers = {"Content-Type": "application/json"}
        return 200, response_headers, json.dumps(response)

    def create_web_acl(self) -> TYPE_RESPONSE:
        """https://docs.aws.amazon.com/waf/latest/APIReference/API_CreateWebACL.html (response syntax section)"""

        scope = self._get_param("Scope")
        if scope == "CLOUDFRONT":
            self.region = GLOBAL_REGION
        name = self._get_param("Name")
        body = json.loads(self.body)
        description = body.get("Description")
        tags = body.get("Tags", [])
        rules = body.get("Rules", [])
        web_acl = self.wafv2_backend.create_web_acl(
            name,
            body["VisibilityConfig"],
            body["DefaultAction"],
            scope,
            description,
            tags,
            rules,
        )
        response = {"Summary": web_acl.to_dict()}
        response_headers = {"Content-Type": "application/json"}
        return 200, response_headers, json.dumps(response)

    def delete_web_acl(self) -> TYPE_RESPONSE:
        scope = self._get_param("Scope")
        if scope == "CLOUDFRONT":
            self.region = GLOBAL_REGION
        name = self._get_param("Name")
        _id = self._get_param("Id")
        self.wafv2_backend.delete_web_acl(name, _id)
        response_headers = {"Content-Type": "application/json"}
        return 200, response_headers, "{}"

    def get_web_acl(self) -> TYPE_RESPONSE:
        scope = self._get_param("Scope")
        if scope == "CLOUDFRONT":
            self.region = GLOBAL_REGION
        name = self._get_param("Name")
        _id = self._get_param("Id")
        web_acl = self.wafv2_backend.get_web_acl(name, _id)
        response = {"WebACL": web_acl.to_dict(), "LockToken": web_acl.lock_token}
        response_headers = {"Content-Type": "application/json"}
        return 200, response_headers, json.dumps(response)

    def list_web_ac_ls(self) -> TYPE_RESPONSE:
        """https://docs.aws.amazon.com/waf/latest/APIReference/API_ListWebACLs.html (response syntax section)"""

        scope = self._get_param("Scope")
        if scope == "CLOUDFRONT":
            self.region = GLOBAL_REGION
        all_web_acls = self.wafv2_backend.list_web_acls()
        response = {"NextMarker": "Not Implemented", "WebACLs": all_web_acls}
        response_headers = {"Content-Type": "application/json"}
        return 200, response_headers, json.dumps(response)

    def list_rule_groups(self) -> TYPE_RESPONSE:
        scope = self._get_param("Scope")
        if scope == "CLOUDFRONT":
            self.region = GLOBAL_REGION
        rule_groups = self.wafv2_backend.list_rule_groups()
        response = {"RuleGroups": [rg.to_dict() for rg in rule_groups]}
        response_headers = {"Content-Type": "application/json"}
        return 200, response_headers, json.dumps(response)

    def list_tags_for_resource(self) -> TYPE_RESPONSE:
        arn = self._get_param("ResourceARN")
        self.region = arn.split(":")[3]
        tags = self.wafv2_backend.list_tags_for_resource(arn)
        response = {"TagInfoForResource": {"ResourceARN": arn, "TagList": tags}}
        response_headers = {"Content-Type": "application/json"}
        return 200, response_headers, json.dumps(response)

    def tag_resource(self) -> TYPE_RESPONSE:
        body = json.loads(self.body)
        arn = body.get("ResourceARN")
        self.region = arn.split(":")[3]
        tags = body.get("Tags")
        self.wafv2_backend.tag_resource(arn, tags)
        return 200, {}, "{}"

    def untag_resource(self) -> TYPE_RESPONSE:
        body = json.loads(self.body)
        arn = body.get("ResourceARN")
        self.region = arn.split(":")[3]
        tag_keys = body.get("TagKeys")
        self.wafv2_backend.untag_resource(arn, tag_keys)
        return 200, {}, "{}"

    def update_web_acl(self) -> TYPE_RESPONSE:
        body = json.loads(self.body)
        name = body.get("Name")
        _id = body.get("Id")
        scope = body.get("Scope")
        if scope == "CLOUDFRONT":
            self.region = GLOBAL_REGION
        default_action = body.get("DefaultAction")
        rules = body.get("Rules")
        description = body.get("Description")
        visibility_config = body.get("VisibilityConfig")
        lock_token = self.wafv2_backend.update_web_acl(
            name, _id, default_action, rules, description, visibility_config
        )
        return 200, {}, json.dumps({"NextLockToken": lock_token})

    def create_ip_set(self) -> TYPE_RESPONSE:
        body = json.loads(self.body)

        name = body.get("Name")
        scope = body.get("Scope")

        if scope == "CLOUDFRONT":
            self.region = GLOBAL_REGION

        description = body.get("Description")
        ip_address_version = body.get("IPAddressVersion")
        addresses = body.get("Addresses")
        tags = body.get("Tags")

        ip_set = self.wafv2_backend.create_ip_set(
            name, scope, description, ip_address_version, addresses, tags
        )
        return (
            200,
            {},
            json.dumps(
                {
                    "Summary": {
                        "Name": ip_set.name,
                        "Id": ip_set.ip_set_id,
                        "Description": ip_set.description,
                        "LockToken": ip_set.lock_token,
                        "ARN": ip_set.arn,
                    }
                }
            ),
        )

    def delete_ip_set(self) -> TYPE_RESPONSE:
        body = json.loads(self.body)

        name = body.get("Name")
        scope = body.get("Scope")
        if scope == "CLOUDFRONT":
            self.region = GLOBAL_REGION

        _id = body.get("Id")
        lock_token = body.get("LockToken")

        self.wafv2_backend.delete_ip_set(name, scope, _id, lock_token)

        return 200, {}, "{}"

    def list_ip_sets(self) -> TYPE_RESPONSE:
        scope = self._get_param("Scope")
        if scope == "CLOUDFRONT":
            self.region = GLOBAL_REGION
        ip_sets = self.wafv2_backend.list_ip_sets(scope)

        formatted_ip_sets = [
            {
                "Name": ip_set.name,
                "Id": ip_set.ip_set_id,
                "Description": ip_set.description,
                "LockToken": ip_set.lock_token,
                "ARN": ip_set.arn,
            }
            for ip_set in ip_sets
        ]

        return (
            200,
            {},
            json.dumps(
                {"NextMarker": str(mock_random.uuid4()), "IPSets": formatted_ip_sets}
            ),
        )

    def get_ip_set(self) -> TYPE_RESPONSE:
        scope = self._get_param("Scope")
        if scope == "CLOUDFRONT":
            self.region = GLOBAL_REGION

        ip_set = self.wafv2_backend.get_ip_set(
            name=self._get_param("Name"), scope=scope, _id=self._get_param("Id")
        )

        dict_ip = ip_set.to_dict()
        lock_token = dict_ip.pop("LockToken")
        return 200, {}, json.dumps({"IPSet": dict_ip, "LockToken": lock_token})

    def update_ip_set(self) -> TYPE_RESPONSE:
        body = json.loads(self.body)

        name = body.get("Name")
        scope = body.get("Scope")
        _id = body.get("Id")

        if scope == "CLOUDFRONT":
            self.region = GLOBAL_REGION

        description = body.get("Description")
        addresses = body.get("Addresses")
        lock_token = body.get("LockToken")
        updated_ip_set = self.wafv2_backend.update_ip_set(
            name, scope, _id, description, addresses, lock_token
        )

        return 200, {}, json.dumps({"NextLockToken": updated_ip_set.lock_token})

    def put_logging_configuration(self) -> TYPE_RESPONSE:
        body = json.loads(self.body)
        logging_configuration_parameter = body["LoggingConfiguration"]
        resource_arn = logging_configuration_parameter["ResourceArn"]
        log_destination_configs = logging_configuration_parameter[
            "LogDestinationConfigs"
        ]
        redacted_fields = logging_configuration_parameter.get("RedactedFields")
        managed_by_firewall_manager = logging_configuration_parameter.get(
            "ManagedByFirewallManager"
        )
        logging_filter = logging_configuration_parameter.get("LoggingFilter")
        logging_configuration = self.wafv2_backend.put_logging_configuration(
            resource_arn,
            log_destination_configs,
            redacted_fields,
            managed_by_firewall_manager,
            logging_filter,
        )
        return (
            200,
            {},
            json.dumps(
                {
                    "LoggingConfiguration": {
                        k: v
                        for k, v in logging_configuration.to_dict().items()
                        if v is not None
                    }
                }
            ),
        )

    def get_logging_configuration(self) -> TYPE_RESPONSE:
        body = json.loads(self.body)
        resource_arn = body["ResourceArn"]
        logging_configuration = self.wafv2_backend.get_logging_configuration(
            resource_arn
        )
        return (
            200,
            {},
            json.dumps(
                {
                    "LoggingConfiguration": {
                        k: v
                        for k, v in logging_configuration.to_dict().items()
                        if v is not None
                    }
                }
            ),
        )

    def list_logging_configurations(self) -> TYPE_RESPONSE:
        body = json.loads(self.body)
        scope = body.get("Scope")
        log_configs = self.wafv2_backend.list_logging_configurations(scope)

        formatted = [
            {k: v for k, v in config.to_dict().items() if v is not None}
            for config in log_configs
        ]
        return 200, {}, json.dumps({"LoggingConfigurations": formatted})

    def delete_logging_configuration(self) -> TYPE_RESPONSE:
        body = json.loads(self.body)
        resource_arn = body["ResourceArn"]
        self.wafv2_backend.delete_logging_configuration(resource_arn)
        return 200, {}, "{}"


# notes about region and scope
# --scope = CLOUDFRONT is ALWAYS us-east-1 (but we use "global" instead to differentiate between REGIONAL us-east-1)
# --scope = REGIONAL defaults to us-east-1, but could be anything if specified with --region=<anyRegion>
# region is grabbed from the auth header, NOT from the body - even with --region flag
# The CLOUDFRONT wacls in aws console are located in us-east-1 but the us-east-1 REGIONAL wacls are not included
