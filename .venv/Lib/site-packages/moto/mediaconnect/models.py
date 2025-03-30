from collections import OrderedDict
from typing import Any, Dict, List, Optional

from moto.core.base_backend import BackendDict, BaseBackend
from moto.core.common_models import BaseModel
from moto.mediaconnect.exceptions import NotFoundException
from moto.moto_api._internal import mock_random as random
from moto.utilities.tagging_service import TaggingService
from moto.utilities.utils import get_partition


class Flow(BaseModel):
    def __init__(self, account_id: str, region_name: str, **kwargs: Any):
        self.id = random.uuid4().hex
        self.availability_zone = kwargs.get("availability_zone")
        self.entitlements = kwargs.get("entitlements", [])
        self.name = kwargs.get("name")
        self.outputs = kwargs.get("outputs", [])
        self.source = kwargs.get("source", {})
        self.source_failover_config = kwargs.get("source_failover_config", {})
        self.sources = kwargs.get("sources", [])
        self.vpc_interfaces = kwargs.get("vpc_interfaces", [])
        self.status: Optional[str] = (
            "STANDBY"  # one of 'STANDBY'|'ACTIVE'|'UPDATING'|'DELETING'|'STARTING'|'STOPPING'|'ERROR'
        )
        self._previous_status: Optional[str] = None
        self.description = "A Moto test flow"
        self.flow_arn = f"arn:{get_partition(region_name)}:mediaconnect:{region_name}:{account_id}:flow:{self.id}:{self.name}"
        self.egress_ip = "127.0.0.1"
        if self.source and not self.sources:
            self.sources = [
                self.source,
            ]

    def to_dict(self, include: Optional[List[str]] = None) -> Dict[str, Any]:
        data = {
            "availabilityZone": self.availability_zone,
            "description": self.description,
            "egressIp": self.egress_ip,
            "entitlements": self.entitlements,
            "flowArn": self.flow_arn,
            "name": self.name,
            "outputs": self.outputs,
            "source": self.source,
            "sourceFailoverConfig": self.source_failover_config,
            "sources": self.sources,
            "status": self.status,
            "vpcInterfaces": self.vpc_interfaces,
        }
        if include:
            new_data = {k: v for k, v in data.items() if k in include}
            if "sourceType" in include:
                new_data["sourceType"] = "OWNED"
            return new_data
        return data

    def resolve_transient_states(self) -> None:
        if self.status in ["STARTING"]:
            self.status = "ACTIVE"
        if self.status in ["STOPPING"]:
            self.status = "STANDBY"
        if self.status in ["UPDATING"]:
            self.status = self._previous_status
            self._previous_status = None


class MediaConnectBackend(BaseBackend):
    def __init__(self, region_name: str, account_id: str):
        super().__init__(region_name, account_id)
        self._flows: Dict[str, Flow] = OrderedDict()
        self.tagger = TaggingService()

    def _add_source_details(
        self,
        source: Optional[Dict[str, Any]],
        flow_id: str,
        ingest_ip: str = "127.0.0.1",
    ) -> None:
        if source:
            source["sourceArn"] = (
                f"arn:{get_partition(self.region_name)}:mediaconnect:{self.region_name}:{self.account_id}:source"
                f":{flow_id}:{source['name']}"
            )
            if not source.get("entitlementArn"):
                source["ingestIp"] = ingest_ip

    def _add_entitlement_details(
        self, entitlement: Optional[Dict[str, Any]], entitlement_id: str
    ) -> None:
        if entitlement:
            entitlement["entitlementArn"] = (
                f"arn:{get_partition(self.region_name)}:mediaconnect:{self.region_name}"
                f":{self.account_id}:entitlement:{entitlement_id}"
                f":{entitlement['name']}"
            )

    def _create_flow_add_details(self, flow: Flow) -> None:
        for index, _source in enumerate(flow.sources):
            self._add_source_details(_source, flow.id, f"127.0.0.{index}")

        for index, output in enumerate(flow.outputs or []):
            if output.get("protocol") in ["srt-listener", "zixi-pull"]:
                output["listenerAddress"] = f"{index}.0.0.0"
            output_id = random.uuid4().hex
            arn = (
                f"arn:{get_partition(self.region_name)}:mediaconnect:{self.region_name}"
                f":{self.account_id}:output:{output_id}:{output['name']}"
            )
            output["outputArn"] = arn

        for _, entitlement in enumerate(flow.entitlements):
            entitlement_id = random.uuid4().hex
            self._add_entitlement_details(entitlement, entitlement_id)

    def create_flow(
        self,
        availability_zone: str,
        entitlements: List[Dict[str, Any]],
        name: str,
        outputs: List[Dict[str, Any]],
        source: Dict[str, Any],
        source_failover_config: Dict[str, Any],
        sources: List[Dict[str, Any]],
        vpc_interfaces: List[Dict[str, Any]],
    ) -> Flow:
        flow = Flow(
            account_id=self.account_id,
            region_name=self.region_name,
            availability_zone=availability_zone,
            entitlements=entitlements,
            name=name,
            outputs=outputs,
            source=source,
            source_failover_config=source_failover_config,
            sources=sources,
            vpc_interfaces=vpc_interfaces,
        )
        self._create_flow_add_details(flow)
        self._flows[flow.flow_arn] = flow
        return flow

    def list_flows(self, max_results: Optional[int]) -> List[Dict[str, Any]]:
        """
        Pagination is not yet implemented
        """
        flows = list(self._flows.values())
        if max_results is not None:
            flows = flows[:max_results]
        return [
            fl.to_dict(
                include=[
                    "availabilityZone",
                    "description",
                    "flowArn",
                    "name",
                    "sourceType",
                    "status",
                ]
            )
            for fl in flows
        ]

    def describe_flow(self, flow_arn: str) -> Flow:
        if flow_arn in self._flows:
            flow = self._flows[flow_arn]
            flow.resolve_transient_states()
            return flow
        raise NotFoundException(message="Flow not found.")

    def delete_flow(self, flow_arn: str) -> Flow:
        if flow_arn in self._flows:
            return self._flows.pop(flow_arn)
        raise NotFoundException(message="Flow not found.")

    def start_flow(self, flow_arn: str) -> Flow:
        if flow_arn in self._flows:
            flow = self._flows[flow_arn]
            flow.status = "STARTING"
            return flow
        raise NotFoundException(message="Flow not found.")

    def stop_flow(self, flow_arn: str) -> Flow:
        if flow_arn in self._flows:
            flow = self._flows[flow_arn]
            flow.status = "STOPPING"
            return flow
        raise NotFoundException(message="Flow not found.")

    def tag_resource(self, resource_arn: str, tags: Dict[str, Any]) -> None:
        tag_list = TaggingService.convert_dict_to_tags_input(tags)
        self.tagger.tag_resource(resource_arn, tag_list)

    def list_tags_for_resource(self, resource_arn: str) -> Dict[str, str]:
        if self.tagger.has_tags(resource_arn):
            return self.tagger.get_tag_dict_for_resource(resource_arn)
        raise NotFoundException(message="Resource not found.")

    def add_flow_vpc_interfaces(
        self, flow_arn: str, vpc_interfaces: List[Dict[str, Any]]
    ) -> Flow:
        if flow_arn in self._flows:
            flow = self._flows[flow_arn]
            flow.vpc_interfaces = vpc_interfaces
            return flow
        raise NotFoundException(message=f"flow with arn={flow_arn} not found")

    def add_flow_outputs(self, flow_arn: str, outputs: List[Dict[str, Any]]) -> Flow:
        if flow_arn in self._flows:
            flow = self._flows[flow_arn]
            flow.outputs = outputs
            return flow
        raise NotFoundException(message=f"flow with arn={flow_arn} not found")

    def remove_flow_vpc_interface(self, flow_arn: str, vpc_interface_name: str) -> None:
        if flow_arn in self._flows:
            flow = self._flows[flow_arn]
            flow.vpc_interfaces = [
                vpc_interface
                for vpc_interface in self._flows[flow_arn].vpc_interfaces
                if vpc_interface["name"] != vpc_interface_name
            ]
        else:
            raise NotFoundException(message=f"flow with arn={flow_arn} not found")

    def remove_flow_output(self, flow_arn: str, output_name: str) -> None:
        if flow_arn in self._flows:
            flow = self._flows[flow_arn]
            flow.outputs = [
                output
                for output in self._flows[flow_arn].outputs
                if output["name"] != output_name
            ]
        else:
            raise NotFoundException(message=f"flow with arn={flow_arn} not found")

    def update_flow_output(
        self,
        flow_arn: str,
        output_arn: str,
        cidr_allow_list: List[str],
        description: str,
        destination: str,
        encryption: Dict[str, str],
        max_latency: int,
        media_stream_output_configuration: List[Dict[str, Any]],
        min_latency: int,
        port: int,
        protocol: str,
        remote_id: str,
        sender_control_port: int,
        sender_ip_address: str,
        smoothing_latency: int,
        stream_id: str,
        vpc_interface_attachment: Dict[str, str],
    ) -> Dict[str, Any]:
        if flow_arn not in self._flows:
            raise NotFoundException(message=f"flow with arn={flow_arn} not found")
        flow = self._flows[flow_arn]
        for output in flow.outputs:
            if output["outputArn"] == output_arn:
                output["cidrAllowList"] = cidr_allow_list
                output["description"] = description
                output["destination"] = destination
                output["encryption"] = encryption
                output["maxLatency"] = max_latency
                output["mediaStreamOutputConfiguration"] = (
                    media_stream_output_configuration
                )
                output["minLatency"] = min_latency
                output["port"] = port
                output["protocol"] = protocol
                output["remoteId"] = remote_id
                output["senderControlPort"] = sender_control_port
                output["senderIpAddress"] = sender_ip_address
                output["smoothingLatency"] = smoothing_latency
                output["streamId"] = stream_id
                output["vpcInterfaceAttachment"] = vpc_interface_attachment
                return output
        raise NotFoundException(message=f"output with arn={output_arn} not found")

    def add_flow_sources(
        self, flow_arn: str, sources: List[Dict[str, Any]]
    ) -> List[Dict[str, Any]]:
        if flow_arn not in self._flows:
            raise NotFoundException(message=f"flow with arn={flow_arn} not found")
        flow = self._flows[flow_arn]
        for source in sources:
            source_id = random.uuid4().hex
            name = source["name"]
            arn = f"arn:{get_partition(self.region_name)}:mediaconnect:{self.region_name}:{self.account_id}:source:{source_id}:{name}"
            source["sourceArn"] = arn
        flow.sources = sources
        return sources

    def update_flow_source(
        self,
        flow_arn: str,
        source_arn: str,
        decryption: str,
        description: str,
        entitlement_arn: str,
        ingest_port: int,
        max_bitrate: int,
        max_latency: int,
        max_sync_buffer: int,
        media_stream_source_configurations: List[Dict[str, Any]],
        min_latency: int,
        protocol: str,
        sender_control_port: int,
        sender_ip_address: str,
        stream_id: str,
        vpc_interface_name: str,
        whitelist_cidr: str,
    ) -> Optional[Dict[str, Any]]:
        if flow_arn not in self._flows:
            raise NotFoundException(message=f"flow with arn={flow_arn} not found")
        flow = self._flows[flow_arn]
        source: Optional[Dict[str, Any]] = next(
            iter(
                [source for source in flow.sources if source["sourceArn"] == source_arn]
            ),
            {},
        )
        if source:
            source["decryption"] = decryption
            source["description"] = description
            source["entitlementArn"] = entitlement_arn
            source["ingestPort"] = ingest_port
            source["maxBitrate"] = max_bitrate
            source["maxLatency"] = max_latency
            source["maxSyncBuffer"] = max_sync_buffer
            source["mediaStreamSourceConfigurations"] = (
                media_stream_source_configurations
            )
            source["minLatency"] = min_latency
            source["protocol"] = protocol
            source["senderControlPort"] = sender_control_port
            source["senderIpAddress"] = sender_ip_address
            source["streamId"] = stream_id
            source["vpcInterfaceName"] = vpc_interface_name
            source["whitelistCidr"] = whitelist_cidr
        return source

    def grant_flow_entitlements(
        self,
        flow_arn: str,
        entitlements: List[Dict[str, Any]],
    ) -> List[Dict[str, Any]]:
        if flow_arn not in self._flows:
            raise NotFoundException(message=f"flow with arn={flow_arn} not found")
        flow = self._flows[flow_arn]
        for entitlement in entitlements:
            entitlement_id = random.uuid4().hex
            name = entitlement["name"]
            arn = f"arn:{get_partition(self.region_name)}:mediaconnect:{self.region_name}:{self.account_id}:entitlement:{entitlement_id}:{name}"
            entitlement["entitlementArn"] = arn

        flow.entitlements += entitlements
        return entitlements

    def revoke_flow_entitlement(self, flow_arn: str, entitlement_arn: str) -> None:
        if flow_arn not in self._flows:
            raise NotFoundException(message=f"flow with arn={flow_arn} not found")
        flow = self._flows[flow_arn]
        for entitlement in flow.entitlements:
            if entitlement_arn == entitlement["entitlementArn"]:
                flow.entitlements.remove(entitlement)
                return
        raise NotFoundException(
            message=f"entitlement with arn={entitlement_arn} not found"
        )

    def update_flow_entitlement(
        self,
        flow_arn: str,
        entitlement_arn: str,
        description: str,
        encryption: Dict[str, str],
        entitlement_status: str,
        name: str,
        subscribers: List[str],
    ) -> Dict[str, Any]:
        if flow_arn not in self._flows:
            raise NotFoundException(message=f"flow with arn={flow_arn} not found")
        flow = self._flows[flow_arn]
        for entitlement in flow.entitlements:
            if entitlement_arn == entitlement["entitlementArn"]:
                entitlement["description"] = description
                entitlement["encryption"] = encryption
                entitlement["entitlementStatus"] = entitlement_status
                entitlement["name"] = name
                entitlement["subscribers"] = subscribers
                return entitlement
        raise NotFoundException(
            message=f"entitlement with arn={entitlement_arn} not found"
        )


mediaconnect_backends = BackendDict(MediaConnectBackend, "mediaconnect")
