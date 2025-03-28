import json
from urllib.parse import unquote

from moto.core.responses import BaseResponse

from .models import MediaConnectBackend, mediaconnect_backends


class MediaConnectResponse(BaseResponse):
    def __init__(self) -> None:
        super().__init__(service_name="mediaconnect")

    @property
    def mediaconnect_backend(self) -> MediaConnectBackend:
        return mediaconnect_backends[self.current_account][self.region]

    def create_flow(self) -> str:
        availability_zone = self._get_param("availabilityZone")
        entitlements = self._get_param("entitlements", [])
        name = self._get_param("name")
        outputs = self._get_param("outputs")
        source = self._get_param("source")
        source_failover_config = self._get_param("sourceFailoverConfig")
        sources = self._get_param("sources")
        vpc_interfaces = self._get_param("vpcInterfaces")
        flow = self.mediaconnect_backend.create_flow(
            availability_zone=availability_zone,
            entitlements=entitlements,
            name=name,
            outputs=outputs,
            source=source,
            source_failover_config=source_failover_config,
            sources=sources,
            vpc_interfaces=vpc_interfaces,
        )
        return json.dumps(dict(flow=flow.to_dict()))

    def list_flows(self) -> str:
        max_results = self._get_int_param("maxResults")
        flows = self.mediaconnect_backend.list_flows(max_results=max_results)
        return json.dumps(dict(flows=flows))

    def describe_flow(self) -> str:
        flow_arn = unquote(self._get_param("flowArn"))
        flow = self.mediaconnect_backend.describe_flow(flow_arn=flow_arn)
        return json.dumps(dict(flow=flow.to_dict()))

    def delete_flow(self) -> str:
        flow_arn = unquote(self._get_param("flowArn"))
        flow = self.mediaconnect_backend.delete_flow(flow_arn=flow_arn)
        return json.dumps(dict(flowArn=flow.flow_arn, status=flow.status))

    def start_flow(self) -> str:
        flow_arn = unquote(self._get_param("flowArn"))
        flow = self.mediaconnect_backend.start_flow(flow_arn=flow_arn)
        return json.dumps(dict(flowArn=flow.flow_arn, status=flow.status))

    def stop_flow(self) -> str:
        flow_arn = unquote(self._get_param("flowArn"))
        flow = self.mediaconnect_backend.stop_flow(flow_arn=flow_arn)
        return json.dumps(dict(flowArn=flow.flow_arn, status=flow.status))

    def tag_resource(self) -> str:
        resource_arn = unquote(self._get_param("resourceArn"))
        tags = self._get_param("tags")
        self.mediaconnect_backend.tag_resource(resource_arn=resource_arn, tags=tags)
        return json.dumps(dict())

    def list_tags_for_resource(self) -> str:
        resource_arn = unquote(self._get_param("resourceArn"))
        tags = self.mediaconnect_backend.list_tags_for_resource(resource_arn)
        return json.dumps(dict(tags=tags))

    def add_flow_vpc_interfaces(self) -> str:
        flow_arn = unquote(self._get_param("flowArn"))
        vpc_interfaces = self._get_param("vpcInterfaces")
        flow = self.mediaconnect_backend.add_flow_vpc_interfaces(
            flow_arn=flow_arn, vpc_interfaces=vpc_interfaces
        )
        return json.dumps(
            dict(flow_arn=flow.flow_arn, vpc_interfaces=flow.vpc_interfaces)
        )

    def remove_flow_vpc_interface(self) -> str:
        flow_arn = unquote(self._get_param("flowArn"))
        vpc_interface_name = unquote(self._get_param("vpcInterfaceName"))
        self.mediaconnect_backend.remove_flow_vpc_interface(
            flow_arn=flow_arn, vpc_interface_name=vpc_interface_name
        )
        return json.dumps(
            dict(flow_arn=flow_arn, vpc_interface_name=vpc_interface_name)
        )

    def add_flow_outputs(self) -> str:
        flow_arn = unquote(self._get_param("flowArn"))
        outputs = self._get_param("outputs")
        flow = self.mediaconnect_backend.add_flow_outputs(
            flow_arn=flow_arn, outputs=outputs
        )
        return json.dumps(dict(flow_arn=flow.flow_arn, outputs=flow.outputs))

    def remove_flow_output(self) -> str:
        flow_arn = unquote(self._get_param("flowArn"))
        output_name = unquote(self._get_param("outputArn"))
        self.mediaconnect_backend.remove_flow_output(
            flow_arn=flow_arn, output_name=output_name
        )
        return json.dumps(dict(flow_arn=flow_arn, output_name=output_name))

    def update_flow_output(self) -> str:
        flow_arn = unquote(self._get_param("flowArn"))
        output_arn = unquote(self._get_param("outputArn"))
        cidr_allow_list = self._get_param("cidrAllowList")
        description = self._get_param("description")
        destination = self._get_param("destination")
        encryption = self._get_param("encryption")
        max_latency = self._get_param("maxLatency")
        media_stream_output_configuration = self._get_param(
            "mediaStreamOutputConfiguration"
        )
        min_latency = self._get_param("minLatency")
        port = self._get_param("port")
        protocol = self._get_param("protocol")
        remote_id = self._get_param("remoteId")
        sender_control_port = self._get_param("senderControlPort")
        sender_ip_address = self._get_param("senderIpAddress")
        smoothing_latency = self._get_param("smoothingLatency")
        stream_id = self._get_param("streamId")
        vpc_interface_attachment = self._get_param("vpcInterfaceAttachment")
        output = self.mediaconnect_backend.update_flow_output(
            flow_arn=flow_arn,
            output_arn=output_arn,
            cidr_allow_list=cidr_allow_list,
            description=description,
            destination=destination,
            encryption=encryption,
            max_latency=max_latency,
            media_stream_output_configuration=media_stream_output_configuration,
            min_latency=min_latency,
            port=port,
            protocol=protocol,
            remote_id=remote_id,
            sender_control_port=sender_control_port,
            sender_ip_address=sender_ip_address,
            smoothing_latency=smoothing_latency,
            stream_id=stream_id,
            vpc_interface_attachment=vpc_interface_attachment,
        )
        return json.dumps(dict(flowArn=flow_arn, output=output))

    def add_flow_sources(self) -> str:
        flow_arn = unquote(self._get_param("flowArn"))
        sources = self._get_param("sources")
        sources = self.mediaconnect_backend.add_flow_sources(
            flow_arn=flow_arn, sources=sources
        )
        return json.dumps(dict(flow_arn=flow_arn, sources=sources))

    def update_flow_source(self) -> str:
        flow_arn = unquote(self._get_param("flowArn"))
        source_arn = unquote(self._get_param("sourceArn"))
        description = self._get_param("description")
        decryption = self._get_param("decryption")
        entitlement_arn = self._get_param("entitlementArn")
        ingest_port = self._get_param("ingestPort")
        max_bitrate = self._get_param("maxBitrate")
        max_latency = self._get_param("maxLatency")
        max_sync_buffer = self._get_param("maxSyncbuffer")
        media_stream_source_configurations = self._get_param(
            "mediaStreamSourceConfigurations"
        )
        min_latency = self._get_param("minLatency")
        protocol = self._get_param("protocol")
        sender_control_port = self._get_param("senderControlPort")
        sender_ip_address = self._get_param("senderIpAddress")
        stream_id = self._get_param("streamId")
        vpc_interface_name = self._get_param("vpcInterfaceName")
        whitelist_cidr = self._get_param("whitelistCidr")
        source = self.mediaconnect_backend.update_flow_source(
            flow_arn=flow_arn,
            source_arn=source_arn,
            decryption=decryption,
            description=description,
            entitlement_arn=entitlement_arn,
            ingest_port=ingest_port,
            max_bitrate=max_bitrate,
            max_latency=max_latency,
            max_sync_buffer=max_sync_buffer,
            media_stream_source_configurations=media_stream_source_configurations,
            min_latency=min_latency,
            protocol=protocol,
            sender_control_port=sender_control_port,
            sender_ip_address=sender_ip_address,
            stream_id=stream_id,
            vpc_interface_name=vpc_interface_name,
            whitelist_cidr=whitelist_cidr,
        )
        return json.dumps(dict(flow_arn=flow_arn, source=source))

    def grant_flow_entitlements(self) -> str:
        flow_arn = unquote(self._get_param("flowArn"))
        entitlements = self._get_param("entitlements")
        entitlements = self.mediaconnect_backend.grant_flow_entitlements(
            flow_arn=flow_arn, entitlements=entitlements
        )
        return json.dumps(dict(flow_arn=flow_arn, entitlements=entitlements))

    def revoke_flow_entitlement(self) -> str:
        flow_arn = unquote(self._get_param("flowArn"))
        entitlement_arn = unquote(self._get_param("entitlementArn"))
        self.mediaconnect_backend.revoke_flow_entitlement(
            flow_arn=flow_arn, entitlement_arn=entitlement_arn
        )
        return json.dumps(dict(flowArn=flow_arn, entitlementArn=entitlement_arn))

    def update_flow_entitlement(self) -> str:
        flow_arn = unquote(self._get_param("flowArn"))
        entitlement_arn = unquote(self._get_param("entitlementArn"))
        description = self._get_param("description")
        encryption = self._get_param("encryption")
        entitlement_status = self._get_param("entitlementStatus")
        name = self._get_param("name")
        subscribers = self._get_param("subscribers")
        entitlement = self.mediaconnect_backend.update_flow_entitlement(
            flow_arn=flow_arn,
            entitlement_arn=entitlement_arn,
            description=description,
            encryption=encryption,
            entitlement_status=entitlement_status,
            name=name,
            subscribers=subscribers,
        )
        return json.dumps(dict(flowArn=flow_arn, entitlement=entitlement))
