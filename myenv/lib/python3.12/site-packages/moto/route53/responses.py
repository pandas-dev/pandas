"""Handles Route53 API requests, invokes method and returns response."""

import re
from typing import Any
from urllib.parse import parse_qs

import xmltodict
from jinja2 import Template

from moto.core.common_types import TYPE_RESPONSE
from moto.core.responses import BaseResponse
from moto.core.utils import iso_8601_datetime_with_milliseconds
from moto.route53.exceptions import InvalidChangeBatch
from moto.route53.models import Route53Backend, route53_backends

XMLNS = "https://route53.amazonaws.com/doc/2013-04-01/"


class Route53(BaseResponse):
    """Handler for Route53 requests and responses."""

    def __init__(self) -> None:
        super().__init__(service_name="route53")

    @staticmethod
    def _convert_to_bool(bool_str: Any) -> bool:  # type: ignore[misc]
        if isinstance(bool_str, bool):
            return bool_str

        if isinstance(bool_str, str):
            return str(bool_str).lower() == "true"

        return False

    @property
    def backend(self) -> Route53Backend:
        return route53_backends[self.current_account][self.partition]

    def create_hosted_zone(self) -> TYPE_RESPONSE:
        vpcid = None
        vpcregion = None
        elements = xmltodict.parse(self.body)
        zone_request = elements["CreateHostedZoneRequest"]
        if "HostedZoneConfig" in zone_request:
            zone_config = zone_request["HostedZoneConfig"]
            comment = zone_config["Comment"]
            if zone_request.get("VPC", {}).get("VPCId", None):
                private_zone = True
            else:
                private_zone = self._convert_to_bool(
                    zone_config.get("PrivateZone", False)
                )
        else:
            comment = None
            private_zone = False

        # It is possible to create a Private Hosted Zone without
        # associating VPC at the time of creation.
        if self._convert_to_bool(private_zone):
            if zone_request.get("VPC", None) is not None:
                vpcid = zone_request["VPC"].get("VPCId", None)
                vpcregion = zone_request["VPC"].get("VPCRegion", None)

        name = zone_request["Name"]
        caller_reference = zone_request["CallerReference"]

        if name[-1] != ".":
            name += "."
        delegation_set_id = zone_request.get("DelegationSetId")

        new_zone = self.backend.create_hosted_zone(
            name,
            comment=comment,
            private_zone=private_zone,
            caller_reference=caller_reference,
            vpcid=vpcid,
            vpcregion=vpcregion,
            delegation_set_id=delegation_set_id,
        )
        template = Template(CREATE_HOSTED_ZONE_RESPONSE).render(zone=new_zone)
        headers = {
            "Location": f"https://route53.amazonaws.com/2013-04-01/hostedzone/{new_zone.id}",
            "status": 201,
        }
        return 201, headers, template

    def list_hosted_zones(self) -> str:
        max_size = self.querystring.get("maxitems", [None])[0]
        if max_size:
            max_size = int(max_size)
        marker = self.querystring.get("marker", [None])[0]
        zone_page, next_marker = self.backend.list_hosted_zones(
            marker=marker, max_size=max_size
        )
        template = Template(LIST_HOSTED_ZONES_RESPONSE).render(
            zones=zone_page,
            marker=marker,
            next_marker=next_marker,
            max_items=max_size,
        )
        return template

    def list_hosted_zones_by_name(self) -> str:
        query_params = parse_qs(self.parsed_url.query)
        dnsnames = query_params.get("dnsname")

        dnsname, zones = self.backend.list_hosted_zones_by_name(dnsnames)

        template = Template(LIST_HOSTED_ZONES_BY_NAME_RESPONSE)
        return template.render(zones=zones, dnsname=dnsname, xmlns=XMLNS)

    def list_hosted_zones_by_vpc(self) -> str:
        query_params = parse_qs(self.parsed_url.query)
        vpc_id = query_params.get("vpcid")[0]  # type: ignore
        zones = self.backend.list_hosted_zones_by_vpc(vpc_id)
        template = Template(LIST_HOSTED_ZONES_BY_VPC_RESPONSE)
        return template.render(zones=zones, xmlns=XMLNS)

    def get_hosted_zone_count(self) -> str:
        num_zones = self.backend.get_hosted_zone_count()
        template = Template(GET_HOSTED_ZONE_COUNT_RESPONSE)
        return template.render(zone_count=num_zones, xmlns=XMLNS)

    def get_hosted_zone(self) -> str:
        zoneid = self.parsed_url.path.rstrip("/").rsplit("/", 1)[1]
        the_zone = self.backend.get_hosted_zone(zoneid)
        template = Template(GET_HOSTED_ZONE_RESPONSE)
        return template.render(zone=the_zone)

    def delete_hosted_zone(self) -> str:
        zoneid = self.parsed_url.path.rstrip("/").rsplit("/", 1)[1]
        self.backend.delete_hosted_zone(zoneid)
        return DELETE_HOSTED_ZONE_RESPONSE

    def update_hosted_zone_comment(self) -> str:
        zoneid = self.parsed_url.path.rstrip("/").rsplit("/", 1)[1]

        elements = xmltodict.parse(self.body)
        comment = elements.get("UpdateHostedZoneCommentRequest", {}).get(
            "Comment", None
        )
        zone = self.backend.update_hosted_zone_comment(zoneid, comment)
        template = Template(UPDATE_HOSTED_ZONE_COMMENT_RESPONSE)
        return template.render(zone=zone)

    def get_dnssec(self) -> str:
        # TODO: implement enable/disable dnssec apis
        zoneid = self.parsed_url.path.rstrip("/").rsplit("/", 2)[1]

        self.backend.get_dnssec(zoneid)
        return GET_DNSSEC

    def associate_vpc_with_hosted_zone(self) -> str:
        zoneid = self.parsed_url.path.rstrip("/").rsplit("/", 2)[1]

        elements = xmltodict.parse(self.body)
        comment = elements.get("AssociateVPCWithHostedZoneRequest", {}).get(
            "Comment", {}
        )
        vpc = elements.get("AssociateVPCWithHostedZoneRequest", {}).get("VPC", {})
        vpcid = vpc.get("VPCId", None)
        vpcregion = vpc.get("VPCRegion", None)

        self.backend.associate_vpc_with_hosted_zone(zoneid, vpcid, vpcregion)

        template = Template(ASSOCIATE_VPC_RESPONSE)
        return template.render(comment=comment)

    def disassociate_vpc_from_hosted_zone(self) -> str:
        zoneid = self.parsed_url.path.rstrip("/").rsplit("/", 2)[1]

        elements = xmltodict.parse(self.body)
        comment = elements.get("DisassociateVPCFromHostedZoneRequest", {}).get(
            "Comment", {}
        )
        vpc = elements.get("DisassociateVPCFromHostedZoneRequest", {}).get("VPC", {})
        vpcid = vpc.get("VPCId", None)

        self.backend.disassociate_vpc_from_hosted_zone(zoneid, vpcid)

        template = Template(DISASSOCIATE_VPC_RESPONSE)
        return template.render(comment=comment)

    def change_resource_record_sets(self) -> str:
        zoneid = self.parsed_url.path.rstrip("/").rsplit("/", 2)[1]

        elements = xmltodict.parse(self.body)

        change_list = elements["ChangeResourceRecordSetsRequest"]["ChangeBatch"][
            "Changes"
        ]["Change"]
        if not isinstance(change_list, list):
            change_list = [
                elements["ChangeResourceRecordSetsRequest"]["ChangeBatch"]["Changes"][
                    "Change"
                ]
            ]

        # Enforce quotas https://docs.aws.amazon.com/Route53/latest/DeveloperGuide/DNSLimitations.html#limits-api-requests-changeresourcerecordsets
        #  - A request cannot contain more than 1,000 ResourceRecord elements. When the value of the Action element is UPSERT, each ResourceRecord element is counted twice.
        effective_rr_count = 0
        for value in change_list:
            record_set = value["ResourceRecordSet"]
            if "ResourceRecords" not in record_set or not record_set["ResourceRecords"]:
                continue
            resource_records = list(record_set["ResourceRecords"].values())[0]
            effective_rr_count += len(resource_records)
            if value["Action"] == "UPSERT":
                effective_rr_count += len(resource_records)
        if effective_rr_count > 1000:
            raise InvalidChangeBatch

        self.backend.change_resource_record_sets(zoneid, change_list)

        return CHANGE_RRSET_RESPONSE

    def list_resource_record_sets(self) -> TYPE_RESPONSE:
        zoneid = self.parsed_url.path.rstrip("/").rsplit("/", 2)[1]
        querystring = parse_qs(self.parsed_url.query)
        template = Template(LIST_RRSET_RESPONSE)
        start_type = querystring.get("type", [None])[0]
        start_name = querystring.get("name", [None])[0]
        max_items = int(querystring.get("maxitems", ["300"])[0])

        if start_type and not start_name:
            return 400, {"status": 400}, "The input is not valid"

        (
            record_sets,
            next_name,
            next_type,
            is_truncated,
        ) = self.backend.list_resource_record_sets(
            zoneid,
            start_type=start_type,  # type: ignore
            start_name=start_name,  # type: ignore
            max_items=max_items,
        )
        r_template = template.render(
            record_sets=record_sets,
            next_name=next_name,
            next_type=next_type,
            max_items=max_items,
            is_truncated=is_truncated,
        )
        return 200, {}, r_template

    def create_health_check(self) -> TYPE_RESPONSE:
        json_body = xmltodict.parse(self.body)["CreateHealthCheckRequest"]
        caller_reference = json_body["CallerReference"]
        config = json_body["HealthCheckConfig"]
        health_check_args = {
            "ip_address": config.get("IPAddress"),
            "port": config.get("Port"),
            "type": config["Type"],
            "resource_path": config.get("ResourcePath"),
            "fqdn": config.get("FullyQualifiedDomainName"),
            "search_string": config.get("SearchString"),
            "request_interval": config.get("RequestInterval"),
            "failure_threshold": config.get("FailureThreshold"),
            "health_threshold": config.get("HealthThreshold"),
            "measure_latency": config.get("MeasureLatency"),
            "inverted": config.get("Inverted"),
            "disabled": config.get("Disabled"),
            "enable_sni": config.get("EnableSNI"),
            "children": config.get("ChildHealthChecks", {}).get("ChildHealthCheck"),
            "regions": config.get("Regions", {}).get("Region"),
        }
        health_check = self.backend.create_health_check(
            caller_reference, health_check_args
        )
        template = Template(CREATE_HEALTH_CHECK_RESPONSE)
        return (
            201,
            {"status": 201},
            template.render(health_check=health_check, xmlns=XMLNS),
        )

    def list_health_checks(self) -> str:
        template = Template(LIST_HEALTH_CHECKS_RESPONSE)
        health_checks = self.backend.list_health_checks()
        return template.render(health_checks=health_checks, xmlns=XMLNS)

    def get_health_check(self) -> str:
        health_check_id = self.parsed_url.path.split("/")[-1]

        health_check = self.backend.get_health_check(health_check_id)
        template = Template(GET_HEALTH_CHECK_RESPONSE)
        return template.render(health_check=health_check)

    def delete_health_check(self) -> str:
        health_check_id = self.parsed_url.path.split("/")[-1]
        self.backend.delete_health_check(health_check_id)
        template = Template(DELETE_HEALTH_CHECK_RESPONSE)
        return template.render(xmlns=XMLNS)

    def update_health_check(self) -> str:
        health_check_id = self.parsed_url.path.split("/")[-1]
        config = xmltodict.parse(self.body)["UpdateHealthCheckRequest"]
        health_check_args = {
            "ip_address": config.get("IPAddress"),
            "port": config.get("Port"),
            "resource_path": config.get("ResourcePath"),
            "fqdn": config.get("FullyQualifiedDomainName"),
            "search_string": config.get("SearchString"),
            "failure_threshold": config.get("FailureThreshold"),
            "health_threshold": config.get("HealthThreshold"),
            "inverted": config.get("Inverted"),
            "disabled": config.get("Disabled"),
            "enable_sni": config.get("EnableSNI"),
            "children": config.get("ChildHealthChecks", {}).get("ChildHealthCheck"),
            "regions": config.get("Regions", {}).get("Region"),
        }
        health_check = self.backend.update_health_check(
            health_check_id, health_check_args
        )
        template = Template(UPDATE_HEALTH_CHECK_RESPONSE)
        return template.render(health_check=health_check)

    def get_health_check_status(self) -> str:
        health_check_match = re.search(
            r"healthcheck/(?P<health_check_id>[^/]+)/status$", self.parsed_url.path
        )
        health_check_id = health_check_match.group("health_check_id")  # type: ignore[union-attr]

        self.backend.get_health_check(health_check_id)
        template = Template(GET_HEALTH_CHECK_STATUS_RESPONSE)
        return template.render(timestamp=iso_8601_datetime_with_milliseconds())

    def not_implemented_response(
        self, request: Any, full_url: str, headers: Any
    ) -> TYPE_RESPONSE:
        self.setup_class(request, full_url, headers)

        action = ""
        if "tags" in full_url:
            action = "tags"
        elif "trafficpolicyinstances" in full_url:
            action = "policies"
        raise NotImplementedError(
            f"The action for {action} has not been implemented for route 53"
        )

    def list_or_change_tags_for_resource_request(  # type: ignore[return]
        self, request: Any, full_url: str, headers: Any
    ) -> TYPE_RESPONSE:
        self.setup_class(request, full_url, headers)

        id_ = self.parsed_url.path.split("/")[-1]
        type_ = self.parsed_url.path.split("/")[-2]

        if request.method == "GET":
            tags = self.backend.list_tags_for_resource(id_)
            template = Template(LIST_TAGS_FOR_RESOURCE_RESPONSE)
            return (
                200,
                headers,
                template.render(resource_type=type_, resource_id=id_, tags=tags),
            )

        if request.method == "POST":
            tags = xmltodict.parse(self.body)["ChangeTagsForResourceRequest"]

            if "AddTags" in tags:
                tags = tags["AddTags"]  # type: ignore
            elif "RemoveTagKeys" in tags:
                tags = tags["RemoveTagKeys"]  # type: ignore

            self.backend.change_tags_for_resource(id_, tags)
            template = Template(CHANGE_TAGS_FOR_RESOURCE_RESPONSE)
            return 200, headers, template.render()

    def get_change(self) -> str:
        change_id = self.parsed_url.path.rstrip("/").rsplit("/", 1)[1]
        template = Template(GET_CHANGE_RESPONSE)
        return template.render(change_id=change_id, xmlns=XMLNS)

    def create_query_logging_config(self) -> TYPE_RESPONSE:
        json_body = xmltodict.parse(self.body)["CreateQueryLoggingConfigRequest"]
        hosted_zone_id = json_body["HostedZoneId"]
        log_group_arn = json_body["CloudWatchLogsLogGroupArn"]

        query_logging_config = self.backend.create_query_logging_config(
            self.region, hosted_zone_id, log_group_arn
        )

        template = Template(CREATE_QUERY_LOGGING_CONFIG_RESPONSE)
        headers = {"Location": query_logging_config.location, "status": 201}
        return (
            201,
            headers,
            template.render(query_logging_config=query_logging_config, xmlns=XMLNS),
        )

    def list_query_logging_configs(self) -> str:
        hosted_zone_id = self._get_param("hostedzoneid")
        next_token = self._get_param("nexttoken")
        max_results = self._get_int_param("maxresults")

        # The paginator picks up named arguments, returns tuple.
        # pylint: disable=unbalanced-tuple-unpacking
        all_configs, next_token = self.backend.list_query_logging_configs(
            hosted_zone_id=hosted_zone_id,
            next_token=next_token,
            max_results=max_results,
        )

        template = Template(LIST_QUERY_LOGGING_CONFIGS_RESPONSE)
        return template.render(
            query_logging_configs=all_configs, next_token=next_token, xmlns=XMLNS
        )

    def get_query_logging_config(self) -> str:
        query_logging_config_id = self.parsed_url.path.rstrip("/").rsplit("/", 1)[1]

        query_logging_config = self.backend.get_query_logging_config(
            query_logging_config_id
        )
        template = Template(GET_QUERY_LOGGING_CONFIG_RESPONSE)
        return template.render(query_logging_config=query_logging_config, xmlns=XMLNS)

    def delete_query_logging_config(self) -> str:
        query_logging_config_id = self.parsed_url.path.rstrip("/").rsplit("/", 1)[1]
        self.backend.delete_query_logging_config(query_logging_config_id)
        return ""

    def list_reusable_delegation_sets(self) -> str:
        delegation_sets = self.backend.list_reusable_delegation_sets()
        template = self.response_template(LIST_REUSABLE_DELEGATION_SETS_TEMPLATE)
        return template.render(
            delegation_sets=delegation_sets,
            marker=None,
            is_truncated=False,
            max_items=100,
        )

    def create_reusable_delegation_set(self) -> TYPE_RESPONSE:
        elements = xmltodict.parse(self.body)
        root_elem = elements["CreateReusableDelegationSetRequest"]
        caller_reference = root_elem.get("CallerReference")
        hosted_zone_id = root_elem.get("HostedZoneId")
        delegation_set = self.backend.create_reusable_delegation_set(
            caller_reference=caller_reference, hosted_zone_id=hosted_zone_id
        )
        template = self.response_template(CREATE_REUSABLE_DELEGATION_SET_TEMPLATE)
        return (
            201,
            {"Location": delegation_set.location, "status": 201},
            template.render(delegation_set=delegation_set),
        )

    def get_reusable_delegation_set(self) -> str:
        ds_id = self.parsed_url.path.rstrip("/").rsplit("/")[-1]
        delegation_set = self.backend.get_reusable_delegation_set(
            delegation_set_id=ds_id
        )
        template = self.response_template(GET_REUSABLE_DELEGATION_SET_TEMPLATE)
        return template.render(delegation_set=delegation_set)

    def delete_reusable_delegation_set(self) -> str:
        ds_id = self.parsed_url.path.rstrip("/").rsplit("/")[-1]
        self.backend.delete_reusable_delegation_set(delegation_set_id=ds_id)
        template = self.response_template(DELETE_REUSABLE_DELEGATION_SET_TEMPLATE)
        return template.render()


LIST_TAGS_FOR_RESOURCE_RESPONSE = """
<ListTagsForResourceResponse xmlns="https://route53.amazonaws.com/doc/2015-01-01/">
    <ResourceTagSet>
        <ResourceType>{{resource_type}}</ResourceType>
        <ResourceId>{{resource_id}}</ResourceId>
        <Tags>
            {% for key, value in tags.items() %}
            <Tag>
                <Key>{{key}}</Key>
                <Value>{{value}}</Value>
            </Tag>
            {% endfor %}
        </Tags>
    </ResourceTagSet>
</ListTagsForResourceResponse>
"""

CHANGE_TAGS_FOR_RESOURCE_RESPONSE = """<ChangeTagsForResourceResponse xmlns="https://route53.amazonaws.com/doc/2015-01-01/">
</ChangeTagsForResourceResponse>
"""

LIST_RRSET_RESPONSE = """<?xml version="1.0" encoding="UTF-8"?>
<ListResourceRecordSetsResponse xmlns="https://route53.amazonaws.com/doc/2012-12-12/">
   <ResourceRecordSets>
       {% for record in record_sets %}
       <ResourceRecordSet>
           <Name>{{ record.name }}</Name>
           <Type>{{ record.type_ }}</Type>
           {% if record.set_identifier %}
               <SetIdentifier>{{ record.set_identifier }}</SetIdentifier>
           {% endif %}
           {% if record.weight %}
               <Weight>{{ record.weight }}</Weight>
           {% endif %}
           {% if record.region %}
               <Region>{{ record.region }}</Region>
           {% endif %}
           {% if record.ttl %}
               <TTL>{{ record.ttl }}</TTL>
           {% endif %}
           {% if record.failover %}
               <Failover>{{ record.failover }}</Failover>
           {% endif %}
           {% if record.multi_value %}
               <MultiValueAnswer>{{ record.multi_value }}</MultiValueAnswer>
           {% endif %}
           {% if record.geo_location %}
           <GeoLocation>
           {% for geo_key in ['ContinentCode','CountryCode','SubdivisionCode'] %}
             {% if record.geo_location[geo_key] %}<{{ geo_key }}>{{ record.geo_location[geo_key] }}</{{ geo_key }}>{% endif %}
           {% endfor %}
           </GeoLocation>
           {% endif %}
           {% if record.alias_target %}
           <AliasTarget>
               <HostedZoneId>{{ record.alias_target['HostedZoneId'] }}</HostedZoneId>
               <DNSName>{{ record.alias_target['DNSName'] }}</DNSName>
               <EvaluateTargetHealth>{{ record.alias_target['EvaluateTargetHealth'] }}</EvaluateTargetHealth>
           </AliasTarget>
           {% else %}
           <ResourceRecords>
               {% for resource in record.records %}
               <ResourceRecord>
                   <Value><![CDATA[{{ resource }}]]></Value>
               </ResourceRecord>
               {% endfor %}
           </ResourceRecords>
           {% endif %}
           {% if record.health_check %}
               <HealthCheckId>{{ record.health_check }}</HealthCheckId>
           {% endif %}
        </ResourceRecordSet>
       {% endfor %}
   </ResourceRecordSets>
   {% if is_truncated %}<NextRecordName>{{ next_name }}</NextRecordName>{% endif %}
   {% if is_truncated %}<NextRecordType>{{ next_type }}</NextRecordType>{% endif %}
   <MaxItems>{{ max_items }}</MaxItems>
   <IsTruncated>{{ 'true' if is_truncated else 'false' }}</IsTruncated>
</ListResourceRecordSetsResponse>"""

CHANGE_RRSET_RESPONSE = """<ChangeResourceRecordSetsResponse xmlns="https://route53.amazonaws.com/doc/2012-12-12/">
   <ChangeInfo>
      <Status>INSYNC</Status>
      <SubmittedAt>2010-09-10T01:36:41.958Z</SubmittedAt>
      <Id>/change/C2682N5HXP0BZ4</Id>
   </ChangeInfo>
</ChangeResourceRecordSetsResponse>"""

DELETE_HOSTED_ZONE_RESPONSE = """<DeleteHostedZoneResponse xmlns="https://route53.amazonaws.com/doc/2012-12-12/">
    <ChangeInfo>
      <Status>INSYNC</Status>
      <SubmittedAt>2010-09-10T01:36:41.958Z</SubmittedAt>
      <Id>/change/C2682N5HXP0BZ4</Id>
   </ChangeInfo>
</DeleteHostedZoneResponse>"""

GET_HOSTED_ZONE_COUNT_RESPONSE = """<GetHostedZoneCountResponse> xmlns="https://route53.amazonaws.com/doc/2012-12-12/">
   <HostedZoneCount>{{ zone_count }}</HostedZoneCount>
</GetHostedZoneCountResponse>"""


GET_HOSTED_ZONE_RESPONSE = """<GetHostedZoneResponse xmlns="https://route53.amazonaws.com/doc/2012-12-12/">
   <HostedZone>
      <Id>/hostedzone/{{ zone.id }}</Id>
      <Name>{{ zone.name }}</Name>
        <CallerReference>{{ zone.caller_reference }}</CallerReference>
      <ResourceRecordSetCount>{{ zone.rrsets|count }}</ResourceRecordSetCount>
      <Config>
        {% if zone.comment %}
            <Comment>{{ zone.comment }}</Comment>
        {% endif %}
        <PrivateZone>{{ 'true' if zone.private_zone else 'false' }}</PrivateZone>
      </Config>
   </HostedZone>
   {% if not zone.private_zone %}
   <DelegationSet>
      <Id>{{ zone.delegation_set.id }}</Id>
      <NameServers>
        {% for name in zone.delegation_set.name_servers %}<NameServer>{{ name }}</NameServer>{% endfor %}
      </NameServers>
   </DelegationSet>
   {% endif %}
   {% if zone.private_zone %}
   <VPCs>
      {% for vpc in zone.vpcs %}
      <VPC>
         <VPCId>{{vpc.vpc_id}}</VPCId>
         <VPCRegion>{{vpc.vpc_region}}</VPCRegion>
      </VPC>
      {% endfor %}
   </VPCs>
   {% endif %}
</GetHostedZoneResponse>"""

CREATE_HOSTED_ZONE_RESPONSE = """<CreateHostedZoneResponse xmlns="https://route53.amazonaws.com/doc/2012-12-12/">
    {% if zone.private_zone %}
    <VPC>
      <VPCId>{{zone.vpcid}}</VPCId>
      <VPCRegion>{{zone.vpcregion}}</VPCRegion>
    </VPC>
    {% endif %}
   <HostedZone>
      <Id>/hostedzone/{{ zone.id }}</Id>
      <Name>{{ zone.name }}</Name>
      <CallerReference>{{ zone.caller_reference }}</CallerReference>
      <ResourceRecordSetCount>{{ zone.rrsets|count }}</ResourceRecordSetCount>
      <Config>
        {% if zone.comment %}
            <Comment>{{ zone.comment }}</Comment>
        {% endif %}
        <PrivateZone>{{ 'true' if zone.private_zone else 'false' }}</PrivateZone>
      </Config>
   </HostedZone>
   {% if not zone.private_zone %}
   <DelegationSet>
      <Id>{{ zone.delegation_set.id }}</Id>
      <NameServers>
         {% for name in zone.delegation_set.name_servers %}<NameServer>{{ name }}</NameServer>{% endfor %}
      </NameServers>
   </DelegationSet>
   {% endif %}
   <ChangeInfo>
      <Id>/change/C1PA6795UKMFR9</Id>
      <Status>INSYNC</Status>
      <SubmittedAt>2017-03-15T01:36:41.958Z</SubmittedAt>
   </ChangeInfo>
</CreateHostedZoneResponse>"""

LIST_HOSTED_ZONES_RESPONSE = """<ListHostedZonesResponse xmlns="https://route53.amazonaws.com/doc/2012-12-12/">
   <HostedZones>
      {% for zone in zones %}
      <HostedZone>
         <Id>/hostedzone/{{ zone.id }}</Id>
         <Name>{{ zone.name }}</Name>
         <CallerReference>{{ zone.caller_reference }}</CallerReference>
         <Config>
            {% if zone.comment %}
                <Comment>{{ zone.comment }}</Comment>
            {% endif %}
           <PrivateZone>{{ 'true' if zone.private_zone else 'false' }}</PrivateZone>
         </Config>
         <ResourceRecordSetCount>{{ zone.rrsets|count  }}</ResourceRecordSetCount>
      </HostedZone>
      {% endfor %}
   </HostedZones>
   {% if marker %}<Marker>{{ marker }}</Marker>{% endif %}
   {%if next_marker %}<NextMarker>{{ next_marker }}</NextMarker>{% endif %}
   {%if max_items %}<MaxItems>{{ max_items }}</MaxItems>{% endif %}
   <IsTruncated>{{ 'true' if next_marker else 'false'}}</IsTruncated>
</ListHostedZonesResponse>"""

LIST_HOSTED_ZONES_BY_NAME_RESPONSE = """<ListHostedZonesByNameResponse xmlns="{{ xmlns }}">
  {% if dnsname %}
  <DNSName>{{ dnsname }}</DNSName>
  {% endif %}
  <HostedZones>
      {% for zone in zones %}
      <HostedZone>
         <Id>/hostedzone/{{ zone.id }}</Id>
         <Name>{{ zone.name }}</Name>
         <CallerReference>{{ zone.caller_reference }}</CallerReference>
         <Config>
            {% if zone.comment %}
                <Comment>{{ zone.comment }}</Comment>
            {% endif %}
           <PrivateZone>{{ 'true' if zone.private_zone else 'false' }}</PrivateZone>
         </Config>
         <ResourceRecordSetCount>{{ zone.rrsets|count  }}</ResourceRecordSetCount>
      </HostedZone>
      {% endfor %}
   </HostedZones>
   <IsTruncated>false</IsTruncated>
</ListHostedZonesByNameResponse>"""

LIST_HOSTED_ZONES_BY_VPC_RESPONSE = """<ListHostedZonesByVpcResponse xmlns="{{xmlns}}">
   <HostedZoneSummaries>
       {% for zone in zones -%}
       <HostedZoneSummary>
           <HostedZoneId>{{zone["HostedZoneId"]}}</HostedZoneId>
           <Name>{{zone["Name"]}}</Name>
           <Owner>
               {% if zone["Owner"]["OwningAccount"] -%}
               <OwningAccount>{{zone["Owner"]["OwningAccount"]}}</OwningAccount>
               {% endif -%}
               {% if zone["Owner"]["OwningService"] -%}
               <OwningService>zone["Owner"]["OwningService"]</OwningService>
               {% endif -%}
           </Owner>
       </HostedZoneSummary>
       {% endfor -%}
   </HostedZoneSummaries>
</ListHostedZonesByVpcResponse>"""

CREATE_HEALTH_CHECK_RESPONSE = """<?xml version="1.0" encoding="UTF-8"?>
<CreateHealthCheckResponse xmlns="{{ xmlns }}">
  {{ health_check.to_xml() }}
</CreateHealthCheckResponse>"""

UPDATE_HEALTH_CHECK_RESPONSE = """<?xml version="1.0" encoding="UTF-8"?>
<UpdateHealthCheckResponse>
  {{ health_check.to_xml() }}
</UpdateHealthCheckResponse>
"""

LIST_HEALTH_CHECKS_RESPONSE = """<?xml version="1.0" encoding="UTF-8"?>
<ListHealthChecksResponse xmlns="{{ xmlns }}">
   <HealthChecks>
   {% for health_check in health_checks %}
      {{ health_check.to_xml() }}
    {% endfor %}
   </HealthChecks>
   <IsTruncated>false</IsTruncated>
   <MaxItems>{{ health_checks|length }}</MaxItems>
</ListHealthChecksResponse>"""

DELETE_HEALTH_CHECK_RESPONSE = """<?xml version="1.0" encoding="UTF-8"?>
    <DeleteHealthCheckResponse xmlns="{{ xmlns }}">
</DeleteHealthCheckResponse>"""

GET_CHANGE_RESPONSE = """<?xml version="1.0" encoding="UTF-8"?>
<GetChangeResponse xmlns="{{ xmlns }}">
   <ChangeInfo>
      <Status>INSYNC</Status>
      <SubmittedAt>2010-09-10T01:36:41.958Z</SubmittedAt>
      <Id>{{ change_id }}</Id>
   </ChangeInfo>
</GetChangeResponse>"""

CREATE_QUERY_LOGGING_CONFIG_RESPONSE = """<?xml version="1.0" encoding="UTF-8"?>
<CreateQueryLoggingConfigResponse xmlns="{{ xmlns }}">
  {{ query_logging_config.to_xml() }}
</CreateQueryLoggingConfigResponse>"""

GET_QUERY_LOGGING_CONFIG_RESPONSE = """<?xml version="1.0" encoding="UTF-8"?>
<CreateQueryLoggingConfigResponse xmlns="{{ xmlns }}">
  {{ query_logging_config.to_xml() }}
</CreateQueryLoggingConfigResponse>"""

LIST_QUERY_LOGGING_CONFIGS_RESPONSE = """<?xml version="1.0" encoding="UTF-8"?>
<ListQueryLoggingConfigsResponse xmlns="{{ xmlns }}">
   <QueryLoggingConfigs>
      {% for query_logging_config in query_logging_configs %}
         {{ query_logging_config.to_xml() }}
      {% endfor %}
   </QueryLoggingConfigs>
   {% if next_token %}
      <NextToken>{{ next_token }}</NextToken>
   {% endif %}
</ListQueryLoggingConfigsResponse>"""


CREATE_REUSABLE_DELEGATION_SET_TEMPLATE = """<CreateReusableDelegationSetResponse>
  <DelegationSet>
      <Id>{{ delegation_set.id }}</Id>
      <CallerReference>{{ delegation_set.caller_reference }}</CallerReference>
      <NameServers>
        {% for name in delegation_set.name_servers %}<NameServer>{{ name }}</NameServer>{% endfor %}
      </NameServers>
  </DelegationSet>
</CreateReusableDelegationSetResponse>
"""


LIST_REUSABLE_DELEGATION_SETS_TEMPLATE = """<ListReusableDelegationSetsResponse>
  <DelegationSets>
    {% for delegation in delegation_sets %}
    <DelegationSet>
  <Id>{{ delegation.id }}</Id>
  <CallerReference>{{ delegation.caller_reference }}</CallerReference>
  <NameServers>
    {% for name in delegation.name_servers %}<NameServer>{{ name }}</NameServer>{% endfor %}
  </NameServers>
</DelegationSet>
    {% endfor %}
  </DelegationSets>
  <Marker>{{ marker }}</Marker>
  <IsTruncated>{{ is_truncated }}</IsTruncated>
  <MaxItems>{{ max_items }}</MaxItems>
</ListReusableDelegationSetsResponse>
"""


DELETE_REUSABLE_DELEGATION_SET_TEMPLATE = """<DeleteReusableDelegationSetResponse>
  <DeleteReusableDelegationSetResponse/>
</DeleteReusableDelegationSetResponse>
"""

GET_REUSABLE_DELEGATION_SET_TEMPLATE = """<GetReusableDelegationSetResponse>
<DelegationSet>
  <Id>{{ delegation_set.id }}</Id>
  <CallerReference>{{ delegation_set.caller_reference }}</CallerReference>
  <NameServers>
    {% for name in delegation_set.name_servers %}<NameServer>{{ name }}</NameServer>{% endfor %}
  </NameServers>
</DelegationSet>
</GetReusableDelegationSetResponse>
"""

GET_DNSSEC = """<?xml version="1.0"?>
<GetDNSSECResponse>
    <Status>
        <ServeSignature>NOT_SIGNING</ServeSignature>
    </Status>
    <KeySigningKeys/>
</GetDNSSECResponse>
"""

GET_HEALTH_CHECK_RESPONSE = """<?xml version="1.0"?>
<GetHealthCheckResponse>
    {{ health_check.to_xml() }}
</GetHealthCheckResponse>
"""

GET_HEALTH_CHECK_STATUS_RESPONSE = """<?xml version="1.0" encoding="UTF-8"?>
<GetHealthCheckStatusResponse>
   <HealthCheckObservations>
      <HealthCheckObservation>
         <IPAddress>127.0.13.37</IPAddress>
         <Region>us-east-1</Region>
         <StatusReport>
            <CheckedTime>{{ timestamp }}</CheckedTime>
            <Status>Success: HTTP Status Code: 200. Resolved IP: 127.0.13.37. OK</Status>
         </StatusReport>
      </HealthCheckObservation>
   </HealthCheckObservations>
</GetHealthCheckStatusResponse>
"""

UPDATE_HOSTED_ZONE_COMMENT_RESPONSE = """<?xml version="1.0" encoding="UTF-8"?>
<UpdateHostedZoneCommentResponse>
   <HostedZone>
      <Config>
         {% if zone.comment %}
         <Comment>{{ zone.comment }}</Comment>
         {% endif %}
         <PrivateZone>{{ 'true' if zone.private_zone else 'false' }}</PrivateZone>
      </Config>
      <Id>/hostedzone/{{ zone.id }}</Id>
      <Name>{{ zone.name }}</Name>
      <ResourceRecordSetCount>{{ zone.rrsets|count }}</ResourceRecordSetCount>
   </HostedZone>
</UpdateHostedZoneCommentResponse>
"""

ASSOCIATE_VPC_RESPONSE = """<?xml version="1.0" encoding="UTF-8"?>
<AssociateVPCWithHostedZoneResponse>
   <ChangeInfo>
      <Comment>{{ comment or "" }}</Comment>
      <Id>/change/a1b2c3d4</Id>
      <Status>INSYNC</Status>
      <SubmittedAt>2017-03-31T01:36:41.958Z</SubmittedAt>
   </ChangeInfo>
</AssociateVPCWithHostedZoneResponse>
"""

DISASSOCIATE_VPC_RESPONSE = """<?xml version="1.0" encoding="UTF-8"?>
<DisassociateVPCFromHostedZoneResponse xmlns="https://route53.amazonaws.com/doc/2013-04-01/">
   <ChangeInfo>
      <Comment>{{ comment or "" }}</Comment>
      <Id>/change/a1b2c3d4</Id>
      <Status>INSYNC</Status>
      <SubmittedAt>2017-03-31T01:36:41.958Z</SubmittedAt>
   </ChangeInfo>
</DisassociateVPCFromHostedZoneResponse>
"""
