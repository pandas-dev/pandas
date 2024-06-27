"""Route53Backend class with methods for supported APIs."""

import copy
import itertools
import re
import string
from collections import defaultdict
from datetime import datetime
from typing import Any, Dict, List, Optional, Tuple

from jinja2 import Template

from moto.core.base_backend import BackendDict, BaseBackend
from moto.core.common_models import BaseModel, CloudFormationModel
from moto.moto_api._internal import mock_random as random
from moto.route53.exceptions import (
    ConflictingDomainExists,
    DnsNameInvalidForZone,
    HostedZoneNotEmpty,
    InvalidActionValue,
    InvalidCloudWatchArn,
    InvalidInput,
    LastVPCAssociation,
    NoSuchCloudWatchLogsLogGroup,
    NoSuchDelegationSet,
    NoSuchHealthCheck,
    NoSuchHostedZone,
    NoSuchQueryLoggingConfig,
    PublicZoneVPCAssociation,
    QueryLoggingConfigAlreadyExists,
    ResourceRecordAlreadyExists,
)
from moto.utilities.paginator import paginate
from moto.utilities.utils import PARTITION_NAMES, get_partition

from .utils import PAGINATION_MODEL

ROUTE53_ID_CHOICE = string.ascii_uppercase + string.digits


def create_route53_zone_id() -> str:
    # New ID's look like this Z1RWWTK7Y8UDDQ
    return "".join([random.choice(ROUTE53_ID_CHOICE) for _ in range(0, 15)])


def create_route53_caller_reference() -> str:
    timestamp = datetime.now().strftime("%H:%M:%S.%f")
    random_string = "".join(random.choice(string.ascii_letters) for _ in range(6))
    return f"{random_string} {timestamp}"


class DelegationSet(BaseModel):
    def __init__(
        self,
        caller_reference: str,
        name_servers: Optional[List[str]],
        delegation_set_id: Optional[str],
    ):
        self.caller_reference = caller_reference
        self.name_servers = name_servers or [
            "ns-2048.awsdns-64.com",
            "ns-2049.awsdns-65.net",
            "ns-2050.awsdns-66.org",
            "ns-2051.awsdns-67.co.uk",
        ]
        self.id = delegation_set_id or "".join(
            [random.choice(ROUTE53_ID_CHOICE) for _ in range(5)]
        )
        self.location = f"https://route53.amazonaws.com/delegationset/{self.id}"


class HealthCheck(CloudFormationModel):
    def __init__(
        self,
        health_check_id: str,
        caller_reference: str,
        health_check_args: Dict[str, Any],
    ):
        self.id = health_check_id
        self.ip_address = health_check_args.get("ip_address")
        self.port = health_check_args.get("port") or 80
        self.type_ = health_check_args.get("type")
        self.resource_path = health_check_args.get("resource_path")
        self.fqdn = health_check_args.get("fqdn")
        self.search_string = health_check_args.get("search_string")
        self.request_interval = health_check_args.get("request_interval") or 30
        self.failure_threshold = health_check_args.get("failure_threshold") or 3
        self.health_threshold = health_check_args.get("health_threshold")
        self.measure_latency = health_check_args.get("measure_latency") or False
        self.inverted = health_check_args.get("inverted") or False
        self.disabled = health_check_args.get("disabled") or False
        self.enable_sni = health_check_args.get("enable_sni") or True
        self.caller_reference = caller_reference
        self.children = None
        self.regions = None

    def set_children(self, children: Any) -> None:
        if children and isinstance(children, list):
            self.children = children
        elif children and isinstance(children, str):
            self.children = [children]  # type: ignore

    def set_regions(self, regions: Any) -> None:
        if regions and isinstance(regions, list):
            self.regions = regions
        elif regions and isinstance(regions, str):
            self.regions = [regions]  # type: ignore

    @property
    def physical_resource_id(self) -> str:
        return self.id

    @staticmethod
    def cloudformation_name_type() -> str:
        return ""

    @staticmethod
    def cloudformation_type() -> str:
        # https://docs.aws.amazon.com/AWSCloudFormation/latest/UserGuide/aws-resource-route53-healthcheck.html
        return "AWS::Route53::HealthCheck"

    @classmethod
    def create_from_cloudformation_json(  # type: ignore[misc]
        cls,
        resource_name: str,
        cloudformation_json: Any,
        account_id: str,
        region_name: str,
        **kwargs: Any,
    ) -> "HealthCheck":
        properties = cloudformation_json["Properties"]["HealthCheckConfig"]
        health_check_args = {
            "ip_address": properties.get("IPAddress"),
            "port": properties.get("Port"),
            "type": properties["Type"],
            "resource_path": properties.get("ResourcePath"),
            "fqdn": properties.get("FullyQualifiedDomainName"),
            "search_string": properties.get("SearchString"),
            "request_interval": properties.get("RequestInterval"),
            "failure_threshold": properties.get("FailureThreshold"),
        }
        backend = route53_backends[account_id][get_partition(region_name)]
        health_check = backend.create_health_check(
            caller_reference=resource_name, health_check_args=health_check_args
        )
        return health_check

    def to_xml(self) -> str:
        template = Template(
            """<HealthCheck>
            <Id>{{ health_check.id }}</Id>
            <CallerReference>{{ health_check.caller_reference }}</CallerReference>
            <HealthCheckConfig>
                {% if health_check.type_ != "CALCULATED" %}
                    <IPAddress>{{ health_check.ip_address }}</IPAddress>
                    <Port>{{ health_check.port }}</Port>
                {% endif %}
                <Type>{{ health_check.type_ }}</Type>
                {% if health_check.resource_path %}
                    <ResourcePath>{{ health_check.resource_path }}</ResourcePath>
                {% endif %}
                {% if health_check.fqdn %}
                    <FullyQualifiedDomainName>{{ health_check.fqdn }}</FullyQualifiedDomainName>
                {% endif %}
                {% if health_check.type_ != "CALCULATED" %}
                    <RequestInterval>{{ health_check.request_interval }}</RequestInterval>
                    <FailureThreshold>{{ health_check.failure_threshold }}</FailureThreshold>
                    <MeasureLatency>{{ health_check.measure_latency }}</MeasureLatency>
                {% endif %}
                {% if health_check.type_ == "CALCULATED" %}
                    <HealthThreshold>{{ health_check.health_threshold }}</HealthThreshold>
                {% endif %}
                <Inverted>{{ health_check.inverted }}</Inverted>
                <Disabled>{{ health_check.disabled }}</Disabled>
                <EnableSNI>{{ health_check.enable_sni }}</EnableSNI>
                {% if health_check.search_string %}
                    <SearchString>{{ health_check.search_string }}</SearchString>
                {% endif %}
                {% if health_check.children %}
                    <ChildHealthChecks>
                    {% for child in health_check.children %}
                        <ChildHealthCheck>{{ child }}</ChildHealthCheck>
                    {% endfor %}
                    </ChildHealthChecks>
                {% endif %}
                {% if health_check.regions %}
                    <Regions>
                    {% for region in health_check.regions %}
                        <Region>{{ region }}</Region>
                    {% endfor %}
                    </Regions>
                {% endif %}
            </HealthCheckConfig>
            <HealthCheckVersion>1</HealthCheckVersion>
        </HealthCheck>"""
        )
        return template.render(health_check=self)


class RecordSet(CloudFormationModel):
    def __init__(self, kwargs: Dict[str, Any]):
        self.name = kwargs.get("Name", "")
        self.type_ = kwargs.get("Type")
        self.ttl = kwargs.get("TTL", 0)
        self.records = kwargs.get("ResourceRecords", [])
        self.set_identifier = kwargs.get("SetIdentifier")
        self.weight = kwargs.get("Weight", 0)
        self.region = kwargs.get("Region")
        self.health_check = kwargs.get("HealthCheckId")
        self.hosted_zone_name = kwargs.get("HostedZoneName")
        self.hosted_zone_id = kwargs.get("HostedZoneId")
        self.alias_target = kwargs.get("AliasTarget", [])
        self.failover = kwargs.get("Failover", [])
        self.geo_location = kwargs.get("GeoLocation", [])
        self.multi_value = kwargs.get("MultiValueAnswer")

    @staticmethod
    def cloudformation_name_type() -> str:
        return "Name"

    @staticmethod
    def cloudformation_type() -> str:
        # https://docs.aws.amazon.com/AWSCloudFormation/latest/UserGuide/aws-resource-route53-recordset.html
        return "AWS::Route53::RecordSet"

    @classmethod
    def create_from_cloudformation_json(  # type: ignore[misc]
        cls,
        resource_name: str,
        cloudformation_json: Any,
        account_id: str,
        region_name: str,
        **kwargs: Any,
    ) -> "RecordSet":
        properties = cloudformation_json["Properties"]

        zone_name = properties.get("HostedZoneName")
        backend = route53_backends[account_id][get_partition(region_name)]
        hosted_zone = backend.get_hosted_zone_by_name(zone_name) if zone_name else None
        if hosted_zone is None:
            hosted_zone = backend.get_hosted_zone(properties["HostedZoneId"])
        record_set = hosted_zone.add_rrset(properties)
        return record_set

    @classmethod
    def update_from_cloudformation_json(  # type: ignore[misc]
        cls,
        original_resource: Any,
        new_resource_name: str,
        cloudformation_json: Any,
        account_id: str,
        region_name: str,
    ) -> "RecordSet":
        cls.delete_from_cloudformation_json(
            original_resource.name, cloudformation_json, account_id, region_name
        )
        return cls.create_from_cloudformation_json(
            new_resource_name, cloudformation_json, account_id, region_name
        )

    @classmethod
    def delete_from_cloudformation_json(  # type: ignore[misc]
        cls,
        resource_name: str,
        cloudformation_json: Any,
        account_id: str,
        region_name: str,
    ) -> None:
        # this will break if you changed the zone the record is in,
        # unfortunately
        properties = cloudformation_json["Properties"]

        zone_name = properties.get("HostedZoneName")
        backend = route53_backends[account_id][get_partition(region_name)]
        hosted_zone = backend.get_hosted_zone_by_name(zone_name) if zone_name else None
        if hosted_zone is None:
            hosted_zone = backend.get_hosted_zone(properties["HostedZoneId"])

        try:
            hosted_zone.delete_rrset({"Name": resource_name})
        except KeyError:
            pass

    @property
    def physical_resource_id(self) -> str:
        return self.name

    def delete(self, account_id: str, region: str) -> None:
        """Not exposed as part of the Route 53 API - used for CloudFormation"""
        backend = route53_backends[account_id][get_partition(region)]
        hosted_zone = (
            backend.get_hosted_zone_by_name(self.hosted_zone_name)
            if self.hosted_zone_name
            else None
        )
        if not hosted_zone:
            hosted_zone = backend.get_hosted_zone(self.hosted_zone_id)  # type: ignore[arg-type]
        hosted_zone.delete_rrset({"Name": self.name, "Type": self.type_})


def reverse_domain_name(domain_name: str) -> str:
    if domain_name.endswith("."):  # normalize without trailing dot
        domain_name = domain_name[:-1]
    return ".".join(reversed(domain_name.split(".")))


class ChangeList(List[Dict[str, Any]]):
    """
    Contains a 'clean' list of ResourceRecordChangeSets
    """

    def append(self, item: Any) -> None:
        item["ResourceRecordSet"]["Name"] = item["ResourceRecordSet"]["Name"].strip(".")
        super().append(item)

    def __contains__(self, item: Any) -> bool:
        item["ResourceRecordSet"]["Name"] = item["ResourceRecordSet"]["Name"].strip(".")
        return super().__contains__(item)

    def has_insert_or_update(self, new_rr_set: Dict[str, Any]) -> bool:
        """
        Check if a CREATE or UPSERT record exists where the name and type is the same as the provided record
        If the existing record has TTL/ResourceRecords, the new TTL should have the same
        """
        for change in self:
            if change["Action"] in ["CREATE", "UPSERT"]:
                rr_set = change["ResourceRecordSet"]
                if (
                    rr_set["Name"] == new_rr_set["Name"].strip(".")
                    and rr_set["Type"] == new_rr_set["Type"]
                ):
                    if "TTL" in rr_set:
                        if rr_set["TTL"] == new_rr_set.get("TTL") and rr_set[
                            "ResourceRecords"
                        ] == new_rr_set.get("ResourceRecords"):
                            return True
                    else:
                        return True
        return False


class FakeZone(CloudFormationModel):
    def __init__(
        self,
        name: str,
        id_: str,
        private_zone: bool,
        caller_reference: str,
        comment: Optional[str] = None,
        delegation_set: Optional[DelegationSet] = None,
    ):
        self.name = name
        self.id = id_
        self.vpcs: List[Dict[str, Any]] = []
        if comment is not None:
            self.comment = comment
        self.caller_reference = caller_reference
        self.private_zone = private_zone
        self.rrsets: List[RecordSet] = []
        self.delegation_set = delegation_set
        self.rr_changes = ChangeList()

    def add_rrset(self, record_set: Dict[str, Any]) -> RecordSet:
        record_set_obj = RecordSet(record_set)
        self.rrsets.append(record_set_obj)
        return record_set_obj

    def upsert_rrset(self, record_set: Dict[str, Any]) -> RecordSet:
        new_rrset = RecordSet(record_set)
        for i, rrset in enumerate(self.rrsets):
            if (
                rrset.name == new_rrset.name
                and rrset.type_ == new_rrset.type_
                and rrset.set_identifier == new_rrset.set_identifier
            ):
                self.rrsets[i] = new_rrset
                break
        else:
            self.rrsets.append(new_rrset)
        return new_rrset

    def delete_rrset(self, rrset: Dict[str, Any]) -> None:
        self.rrsets = [
            record_set
            for record_set in self.rrsets
            if record_set.name != rrset["Name"]
            or (rrset.get("Type") is not None and record_set.type_ != rrset["Type"])
        ]

    def delete_rrset_by_id(self, set_identifier: str) -> None:
        self.rrsets = [
            record_set
            for record_set in self.rrsets
            if record_set.set_identifier != set_identifier
        ]

    def add_vpc(
        self, vpc_id: Optional[str], vpc_region: Optional[str]
    ) -> Dict[str, Any]:
        vpc = {}
        if vpc_id is not None:
            vpc["vpc_id"] = vpc_id
        if vpc_region is not None:
            vpc["vpc_region"] = vpc_region
        if vpc_id or vpc_region:
            self.vpcs.append(vpc)
        return vpc

    def delete_vpc(self, vpc_id: str) -> None:
        self.vpcs = [vpc for vpc in self.vpcs if vpc["vpc_id"] != vpc_id]

    def get_record_sets(self, start_type: str, start_name: str) -> List[RecordSet]:
        def predicate(rrset: RecordSet) -> bool:
            rrset_name_reversed = reverse_domain_name(rrset.name)
            start_name_reversed = reverse_domain_name(start_name)
            return rrset_name_reversed < start_name_reversed or (
                rrset_name_reversed == start_name_reversed and rrset.type_ < start_type  # type: ignore
            )

        record_sets = sorted(
            self.rrsets,
            key=lambda rrset: (reverse_domain_name(rrset.name), rrset.type_),
        )

        if start_name:
            start_type = start_type or ""
            record_sets = itertools.dropwhile(predicate, record_sets)  # type: ignore

        return record_sets

    @property
    def physical_resource_id(self) -> str:
        return self.id

    @staticmethod
    def cloudformation_name_type() -> str:
        return "Name"

    @staticmethod
    def cloudformation_type() -> str:
        # https://docs.aws.amazon.com/AWSCloudFormation/latest/UserGuide/aws-resource-route53-hostedzone.html
        return "AWS::Route53::HostedZone"

    @classmethod
    def create_from_cloudformation_json(  # type: ignore[misc]
        cls,
        resource_name: str,
        cloudformation_json: Any,
        account_id: str,
        region_name: str,
        **kwargs: Any,
    ) -> "FakeZone":
        hosted_zone = route53_backends[account_id][
            get_partition(region_name)
        ].create_hosted_zone(resource_name, private_zone=False)
        return hosted_zone


class RecordSetGroup(CloudFormationModel):
    def __init__(self, region_name: str, hosted_zone_id: str, record_sets: List[str]):
        self.region_name = region_name
        self.hosted_zone_id = hosted_zone_id
        self.record_sets = record_sets

    @property
    def physical_resource_id(self) -> str:
        return f"arn:{get_partition(self.region_name)}:route53:::hostedzone/{self.hosted_zone_id}"

    @staticmethod
    def cloudformation_name_type() -> str:
        return ""

    @staticmethod
    def cloudformation_type() -> str:
        # https://docs.aws.amazon.com/AWSCloudFormation/latest/UserGuide/aws-resource-route53-recordsetgroup.html
        return "AWS::Route53::RecordSetGroup"

    @classmethod
    def create_from_cloudformation_json(  # type: ignore[misc]
        cls,
        resource_name: str,
        cloudformation_json: Any,
        account_id: str,
        region_name: str,
        **kwargs: Any,
    ) -> "RecordSetGroup":
        properties = cloudformation_json["Properties"]

        zone_name = properties.get("HostedZoneName")
        backend = route53_backends[account_id][get_partition(region_name)]
        hosted_zone = backend.get_hosted_zone_by_name(zone_name) if zone_name else None
        if hosted_zone is None:
            hosted_zone = backend.get_hosted_zone(properties["HostedZoneId"])
        record_sets = properties["RecordSets"]
        for record_set in record_sets:
            hosted_zone.add_rrset(record_set)

        return RecordSetGroup(region_name, hosted_zone.id, record_sets)


class QueryLoggingConfig(BaseModel):
    """QueryLoggingConfig class; this object isn't part of Cloudformation."""

    def __init__(
        self,
        query_logging_config_id: str,
        hosted_zone_id: str,
        cloudwatch_logs_log_group_arn: str,
    ):
        self.hosted_zone_id = hosted_zone_id
        self.cloudwatch_logs_log_group_arn = cloudwatch_logs_log_group_arn
        self.query_logging_config_id = query_logging_config_id
        self.location = f"https://route53.amazonaws.com/2013-04-01/queryloggingconfig/{self.query_logging_config_id}"

    def to_xml(self) -> str:
        template = Template(
            """<QueryLoggingConfig>
                <CloudWatchLogsLogGroupArn>{{ query_logging_config.cloudwatch_logs_log_group_arn }}</CloudWatchLogsLogGroupArn>
                <HostedZoneId>{{ query_logging_config.hosted_zone_id }}</HostedZoneId>
                <Id>{{ query_logging_config.query_logging_config_id }}</Id>
            </QueryLoggingConfig>"""
        )
        # The "Location" value must be put into the header; that's done in
        # responses.py.
        return template.render(query_logging_config=self)


class Route53Backend(BaseBackend):
    def __init__(self, region_name: str, account_id: str):
        super().__init__(region_name, account_id)
        self.zones: Dict[str, FakeZone] = {}
        self.health_checks: Dict[str, HealthCheck] = {}
        self.resource_tags: Dict[str, Any] = defaultdict(dict)
        self.query_logging_configs: Dict[str, QueryLoggingConfig] = {}
        self.delegation_sets: Dict[str, DelegationSet] = dict()

    def _has_prev_conflicting_domain(
        self, name: str, delegation_set_id: Optional[str]
    ) -> bool:
        """Check if a conflicting domain exists in the backend"""
        if not delegation_set_id:
            return False
        for zone in self.zones.values():
            if not zone.delegation_set or zone.delegation_set.id != delegation_set_id:
                # Delegation sets don't match, so can't possibly conflict
                continue
            if (
                zone.name == name
                or zone.name.endswith(f".{name}")
                or name.endswith(f".{zone.name}")
            ):
                return True
        return False

    def create_hosted_zone(
        self,
        name: str,
        private_zone: bool,
        caller_reference: Optional[str] = None,
        vpcid: Optional[str] = None,
        vpcregion: Optional[str] = None,
        comment: Optional[str] = None,
        delegation_set_id: Optional[str] = None,
    ) -> FakeZone:
        if self._has_prev_conflicting_domain(name, delegation_set_id):
            raise ConflictingDomainExists(name, delegation_set_id)
        new_id = create_route53_zone_id()
        caller_reference = caller_reference or create_route53_caller_reference()
        delegation_set = self.create_reusable_delegation_set(
            caller_reference=f"DelSet_{name}", delegation_set_id=delegation_set_id
        )
        # default delegation set does not contains id
        if not delegation_set_id:
            delegation_set.id = ""
        new_zone = FakeZone(
            name,
            new_id,
            caller_reference=caller_reference,
            private_zone=private_zone,
            comment=comment,
            delegation_set=delegation_set,
        )
        # For each public hosted zone that you create, Amazon Route 53 automatically creates a name server (NS) record
        # and a start of authority (SOA) record.
        # https://docs.aws.amazon.com/Route53/latest/DeveloperGuide/SOA-NSrecords.html
        soa_record_set = {
            "Name": f"{name}" + ("" if name.endswith(".") else "."),
            "Type": "SOA",
            "TTL": 900,
            "ResourceRecords": [
                {
                    "Value": f"{delegation_set.name_servers[0]}. hostmaster.example.com. 1 7200 900 1209600 86400"
                }
            ],
        }
        # default nameservers are also part of rrset
        ns_record_set = {
            "Name": name,
            "ResourceRecords": delegation_set.name_servers,
            "TTL": "172800",
            "Type": "NS",
        }
        new_zone.add_rrset(ns_record_set)
        new_zone.add_rrset(soa_record_set)
        new_zone.add_vpc(vpcid, vpcregion)
        self.zones[new_id] = new_zone
        return new_zone

    def get_dnssec(self, zone_id: str) -> None:
        # check if hosted zone exists
        self.get_hosted_zone(zone_id)

    def associate_vpc_with_hosted_zone(
        self, zone_id: str, vpcid: str, vpcregion: str
    ) -> FakeZone:
        zone = self.get_hosted_zone(zone_id)
        if not zone.private_zone:
            raise PublicZoneVPCAssociation()
        zone.add_vpc(vpcid, vpcregion)
        return zone

    def disassociate_vpc_from_hosted_zone(self, zone_id: str, vpcid: str) -> FakeZone:
        zone = self.get_hosted_zone(zone_id)
        if len(zone.vpcs) <= 1:
            raise LastVPCAssociation()
        zone.delete_vpc(vpcid)
        return zone

    def change_tags_for_resource(self, resource_id: str, tags: Any) -> None:
        if "Tag" in tags:
            if isinstance(tags["Tag"], list):
                for tag in tags["Tag"]:
                    self.resource_tags[resource_id][tag["Key"]] = tag["Value"]
            else:
                key, value = (tags["Tag"]["Key"], tags["Tag"]["Value"])
                self.resource_tags[resource_id][key] = value
        else:
            if "Key" in tags:
                if isinstance(tags["Key"], list):
                    for key in tags["Key"]:
                        del self.resource_tags[resource_id][key]
                else:
                    del self.resource_tags[resource_id][tags["Key"]]

    def list_tags_for_resource(self, resource_id: str) -> Dict[str, str]:
        if resource_id in self.resource_tags:
            return self.resource_tags[resource_id]
        return {}

    def list_resource_record_sets(
        self, zone_id: str, start_type: str, start_name: str, max_items: int
    ) -> Tuple[List[RecordSet], Optional[str], Optional[str], bool]:
        """
        The StartRecordIdentifier-parameter is not yet implemented
        """
        the_zone = self.get_hosted_zone(zone_id)
        all_records = list(the_zone.get_record_sets(start_type, start_name))
        records = all_records[0:max_items]
        next_record = all_records[max_items] if len(all_records) > max_items else None
        next_start_name = next_record.name if next_record else None
        next_start_type = next_record.type_ if next_record else None
        is_truncated = next_record is not None
        return records, next_start_name, next_start_type, is_truncated

    def change_resource_record_sets(
        self, zoneid: str, change_list: List[Dict[str, Any]]
    ) -> None:
        the_zone = self.get_hosted_zone(zoneid)

        for value in change_list:
            if value["Action"] == "CREATE" and value in the_zone.rr_changes:
                name = value["ResourceRecordSet"]["Name"] + "."
                _type = value["ResourceRecordSet"]["Type"]
                raise ResourceRecordAlreadyExists(name=name, _type=_type)

        for value in change_list:
            if value["Action"] == "DELETE":
                # To delete a resource record set, you must specify all the same values that you specified when you created it.
                if not the_zone.rr_changes.has_insert_or_update(
                    value["ResourceRecordSet"]
                ):
                    msg = f"Invalid request: Expected exactly one of [AliasTarget, all of [TTL, and ResourceRecords], or TrafficPolicyInstanceId], but found none in Change with [Action=DELETE, Name={value['ResourceRecordSet']['Name']}, Type={value['ResourceRecordSet']['Type']}, SetIdentifier={value['ResourceRecordSet'].get('SetIdentifier', 'null')}]"
                    raise InvalidInput(msg)

        for value in change_list:
            original_change = copy.deepcopy(value)
            action = value["Action"]

            if action not in ("CREATE", "UPSERT", "DELETE"):
                raise InvalidActionValue(action)

            record_set = value["ResourceRecordSet"]

            cleaned_record_name = record_set["Name"].strip(".")
            cleaned_hosted_zone_name = the_zone.name.strip(".")

            if not cleaned_record_name.endswith(cleaned_hosted_zone_name):
                raise DnsNameInvalidForZone(
                    name=record_set["Name"], zone_name=the_zone.name
                )

            if not record_set["Name"].endswith("."):
                record_set["Name"] += "."

            if action in ("CREATE", "UPSERT"):
                if "ResourceRecords" in record_set:
                    resource_records = list(record_set["ResourceRecords"].values())[0]
                    if not isinstance(resource_records, list):
                        # Depending on how many records there are, this may
                        # or may not be a list
                        resource_records = [resource_records]
                    record_set["ResourceRecords"] = [
                        x["Value"] for x in resource_records
                    ]
                if action == "CREATE":
                    the_zone.add_rrset(record_set)
                else:
                    the_zone.upsert_rrset(record_set)
            elif action == "DELETE":
                if "SetIdentifier" in record_set:
                    the_zone.delete_rrset_by_id(record_set["SetIdentifier"])
                else:
                    the_zone.delete_rrset(record_set)
            the_zone.rr_changes.append(original_change)

    @paginate(pagination_model=PAGINATION_MODEL)
    def list_hosted_zones(self) -> List[FakeZone]:
        """
        The parameters DelegationSetId and HostedZoneType are not yet implemented
        """
        return list(self.zones.values())

    def list_hosted_zones_by_name(
        self, dnsnames: Optional[List[str]]
    ) -> Tuple[Optional[str], List[FakeZone]]:
        if dnsnames:
            dnsname = dnsnames[0]
            if dnsname[-1] != ".":
                dnsname += "."
            zones = [zone for zone in self.zones.values() if zone.name == dnsname]
        else:
            dnsname = None
            # sort by names, but with domain components reversed
            # see http://boto3.readthedocs.io/en/latest/reference/services/route53.html#Route53.Client.list_hosted_zones_by_name

            def sort_key(zone: FakeZone) -> str:
                domains = zone.name.split(".")
                if domains[-1] == "":
                    domains = domains[-1:] + domains[:-1]
                return ".".join(reversed(domains))

            zones = sorted(self.zones.values(), key=sort_key)
        return dnsname, zones

    def list_hosted_zones_by_vpc(self, vpc_id: str) -> List[Dict[str, Any]]:
        """
        Pagination is not yet implemented
        """
        zone_list = []
        for zone in self.zones.values():
            if zone.private_zone is True:
                this_zone = self.get_hosted_zone(zone.id)
                for vpc in this_zone.vpcs:
                    if vpc["vpc_id"] == vpc_id:
                        zone_list.append(
                            {
                                "HostedZoneId": zone.id,
                                "Name": zone.name,
                                "Owner": {"OwningAccount": self.account_id},
                            }
                        )

        return zone_list

    def get_hosted_zone(self, id_: str) -> FakeZone:
        the_zone = self.zones.get(id_.replace("/hostedzone/", ""))
        if not the_zone:
            raise NoSuchHostedZone(id_)
        return the_zone

    def get_hosted_zone_count(self) -> int:
        return len(self.zones.values())

    def get_hosted_zone_by_name(self, name: str) -> Optional[FakeZone]:
        for zone in self.zones.values():
            if zone.name == name:
                return zone
        return None

    def delete_hosted_zone(self, id_: str) -> Optional[FakeZone]:
        # Verify it exists
        zone = self.get_hosted_zone(id_)
        if len(zone.rrsets) > 0:
            for rrset in zone.rrsets:
                if rrset.type_ != "NS" and rrset.type_ != "SOA":
                    raise HostedZoneNotEmpty()
        return self.zones.pop(id_.replace("/hostedzone/", ""), None)

    def update_hosted_zone_comment(self, id_: str, comment: str) -> FakeZone:
        zone = self.get_hosted_zone(id_)
        zone.comment = comment
        return zone

    def create_health_check(
        self, caller_reference: str, health_check_args: Dict[str, Any]
    ) -> HealthCheck:
        health_check_id = str(random.uuid4())
        health_check = HealthCheck(health_check_id, caller_reference, health_check_args)
        health_check.set_children(health_check_args.get("children"))
        health_check.set_regions(health_check_args.get("regions"))
        self.health_checks[health_check_id] = health_check
        return health_check

    def update_health_check(
        self, health_check_id: str, health_check_args: Dict[str, Any]
    ) -> HealthCheck:
        health_check = self.health_checks.get(health_check_id)
        if not health_check:
            raise NoSuchHealthCheck(health_check_id)

        if health_check_args.get("ip_address"):
            health_check.ip_address = health_check_args.get("ip_address")
        if health_check_args.get("port"):
            health_check.port = health_check_args.get("port")
        if health_check_args.get("resource_path"):
            health_check.resource_path = health_check_args.get("resource_path")
        if health_check_args.get("fqdn"):
            health_check.fqdn = health_check_args.get("fqdn")
        if health_check_args.get("search_string"):
            health_check.search_string = health_check_args.get("search_string")
        if health_check_args.get("request_interval"):
            health_check.request_interval = health_check_args.get("request_interval")
        if health_check_args.get("failure_threshold"):
            health_check.failure_threshold = health_check_args.get("failure_threshold")
        if health_check_args.get("health_threshold"):
            health_check.health_threshold = health_check_args.get("health_threshold")
        if health_check_args.get("inverted"):
            health_check.inverted = health_check_args.get("inverted")
        if health_check_args.get("disabled"):
            health_check.disabled = health_check_args.get("disabled")
        if health_check_args.get("enable_sni"):
            health_check.enable_sni = health_check_args.get("enable_sni")
        if health_check_args.get("children"):
            health_check.set_children(health_check_args.get("children"))
        if health_check_args.get("regions"):
            health_check.set_regions(health_check_args.get("regions"))

        return health_check

    def list_health_checks(self) -> List[HealthCheck]:
        return list(self.health_checks.values())

    def delete_health_check(self, health_check_id: str) -> None:
        self.health_checks.pop(health_check_id, None)

    def get_health_check(self, health_check_id: str) -> HealthCheck:
        health_check = self.health_checks.get(health_check_id)
        if not health_check:
            raise NoSuchHealthCheck(health_check_id)
        return health_check

    def get_health_check_status(self) -> None:
        pass  # Logic implemented in responses.py

    @staticmethod
    def _validate_arn(region: str, arn: str) -> None:
        match = re.match(
            rf"arn:{get_partition(region)}:logs:{region}:\d{{12}}:log-group:.+", arn
        )
        if not arn or not match:
            raise InvalidCloudWatchArn()

        # The CloudWatch Logs log group must be in the "us-east-1" region.
        match = re.match(r"^(?:[^:]+:){3}(?P<region>[^:]+).*", arn)
        if not match or match.group("region") != "us-east-1":
            raise InvalidCloudWatchArn()

    def create_query_logging_config(
        self, region: str, hosted_zone_id: str, log_group_arn: str
    ) -> QueryLoggingConfig:
        """Process the create_query_logging_config request."""
        # Does the hosted_zone_id exist?
        zones = list(self.zones.values())
        for zone in zones:
            if zone.id == hosted_zone_id:
                break
        else:
            raise NoSuchHostedZone(hosted_zone_id)

        # Ensure CloudWatch Logs log ARN is valid, otherwise raise an error.
        self._validate_arn(region, log_group_arn)

        # Note:  boto3 checks the resource policy permissions before checking
        # whether the log group exists.  moto doesn't have a way of checking
        # the resource policy, so in some instances moto will complain
        # about a log group that doesn't exist whereas boto3 will complain
        # that "The resource policy that you're using for Route 53 query
        # logging doesn't grant Route 53 sufficient permission to create
        # a log stream in the specified log group."

        from moto.logs import logs_backends  # pylint: disable=import-outside-toplevel

        log_groups = logs_backends[self.account_id][region].describe_log_groups()
        for entry in log_groups[0] if log_groups else []:
            if log_group_arn == entry["arn"]:
                break
        else:
            # There is no CloudWatch Logs log group with the specified ARN.
            raise NoSuchCloudWatchLogsLogGroup()

        # Verify there is no existing query log config using the same hosted
        # zone.
        for query_log in self.query_logging_configs.values():
            if query_log.hosted_zone_id == hosted_zone_id:
                raise QueryLoggingConfigAlreadyExists()

        # Create an instance of the query logging config.
        query_logging_config_id = str(random.uuid4())
        query_logging_config = QueryLoggingConfig(
            query_logging_config_id, hosted_zone_id, log_group_arn
        )
        self.query_logging_configs[query_logging_config_id] = query_logging_config
        return query_logging_config

    def delete_query_logging_config(self, query_logging_config_id: str) -> None:
        """Delete query logging config, if it exists."""
        if query_logging_config_id not in self.query_logging_configs:
            raise NoSuchQueryLoggingConfig()
        self.query_logging_configs.pop(query_logging_config_id)

    def get_query_logging_config(
        self, query_logging_config_id: str
    ) -> QueryLoggingConfig:
        """Return query logging config, if it exists."""
        if query_logging_config_id not in self.query_logging_configs:
            raise NoSuchQueryLoggingConfig()
        return self.query_logging_configs[query_logging_config_id]

    @paginate(pagination_model=PAGINATION_MODEL)
    def list_query_logging_configs(
        self, hosted_zone_id: Optional[str] = None
    ) -> List[QueryLoggingConfig]:
        """Return a list of query logging configs."""
        if hosted_zone_id:
            # Does the hosted_zone_id exist?
            zones = list(self.zones.values())
            for zone in zones:
                if zone.id == hosted_zone_id:
                    break
            else:
                raise NoSuchHostedZone(hosted_zone_id)

        return list(self.query_logging_configs.values())

    def create_reusable_delegation_set(
        self,
        caller_reference: str,
        delegation_set_id: Optional[str] = None,
        hosted_zone_id: Optional[str] = None,
    ) -> DelegationSet:
        name_servers: Optional[List[str]] = None
        if hosted_zone_id:
            hosted_zone = self.get_hosted_zone(hosted_zone_id)
            name_servers = hosted_zone.delegation_set.name_servers  # type: ignore
        delegation_set = DelegationSet(
            caller_reference, name_servers, delegation_set_id
        )
        self.delegation_sets[delegation_set.id] = delegation_set
        return delegation_set

    def list_reusable_delegation_sets(self) -> List[DelegationSet]:
        """
        Pagination is not yet implemented
        """
        return list(self.delegation_sets.values())

    def delete_reusable_delegation_set(self, delegation_set_id: str) -> None:
        self.delegation_sets.pop(delegation_set_id, None)

    def get_reusable_delegation_set(self, delegation_set_id: str) -> DelegationSet:
        if delegation_set_id not in self.delegation_sets:
            raise NoSuchDelegationSet(delegation_set_id)
        return self.delegation_sets[delegation_set_id]


route53_backends = BackendDict(
    Route53Backend,
    "route53",
    use_boto3_regions=False,
    additional_regions=PARTITION_NAMES,
)
