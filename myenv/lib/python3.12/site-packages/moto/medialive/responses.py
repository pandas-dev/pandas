import json

from moto.core.responses import BaseResponse

from .models import MediaLiveBackend, medialive_backends


class MediaLiveResponse(BaseResponse):
    def __init__(self) -> None:
        super().__init__(service_name="medialive")

    @property
    def medialive_backend(self) -> MediaLiveBackend:
        return medialive_backends[self.current_account][self.region]

    def create_channel(self) -> str:
        cdi_input_specification = self._get_param("cdiInputSpecification")
        channel_class = self._get_param("channelClass")
        destinations = self._get_param("destinations")
        encoder_settings = self._get_param("encoderSettings")
        input_attachments = self._get_param("inputAttachments")
        input_specification = self._get_param("inputSpecification")
        log_level = self._get_param("logLevel")
        name = self._get_param("name")
        role_arn = self._get_param("roleArn")
        tags = self._get_param("tags")
        channel = self.medialive_backend.create_channel(
            cdi_input_specification=cdi_input_specification,
            channel_class=channel_class,
            destinations=destinations,
            encoder_settings=encoder_settings,
            input_attachments=input_attachments,
            input_specification=input_specification,
            log_level=log_level,
            name=name,
            role_arn=role_arn,
            tags=tags,
        )

        return json.dumps(
            dict(channel=channel.to_dict(exclude=["pipelinesRunningCount"]))
        )

    def list_channels(self) -> str:
        max_results = self._get_int_param("maxResults")
        channels = self.medialive_backend.list_channels(max_results=max_results)

        return json.dumps(dict(channels=channels, nextToken=None))

    def describe_channel(self) -> str:
        channel_id = self._get_param("channelId")
        channel = self.medialive_backend.describe_channel(channel_id=channel_id)
        return json.dumps(channel.to_dict())

    def delete_channel(self) -> str:
        channel_id = self._get_param("channelId")
        channel = self.medialive_backend.delete_channel(channel_id=channel_id)
        return json.dumps(channel.to_dict())

    def start_channel(self) -> str:
        channel_id = self._get_param("channelId")
        channel = self.medialive_backend.start_channel(channel_id=channel_id)
        return json.dumps(channel.to_dict())

    def stop_channel(self) -> str:
        channel_id = self._get_param("channelId")
        channel = self.medialive_backend.stop_channel(channel_id=channel_id)
        return json.dumps(channel.to_dict())

    def update_channel(self) -> str:
        channel_id = self._get_param("channelId")
        cdi_input_specification = self._get_param("cdiInputSpecification")
        destinations = self._get_param("destinations")
        encoder_settings = self._get_param("encoderSettings")
        input_attachments = self._get_param("inputAttachments")
        input_specification = self._get_param("inputSpecification")
        log_level = self._get_param("logLevel")
        name = self._get_param("name")
        role_arn = self._get_param("roleArn")
        channel = self.medialive_backend.update_channel(
            channel_id=channel_id,
            cdi_input_specification=cdi_input_specification,
            destinations=destinations,
            encoder_settings=encoder_settings,
            input_attachments=input_attachments,
            input_specification=input_specification,
            log_level=log_level,
            name=name,
            role_arn=role_arn,
        )
        return json.dumps(dict(channel=channel.to_dict()))

    def create_input(self) -> str:
        destinations = self._get_param("destinations")
        input_devices = self._get_param("inputDevices")
        input_security_groups = self._get_param("inputSecurityGroups")
        media_connect_flows = self._get_param("mediaConnectFlows")
        name = self._get_param("name")
        role_arn = self._get_param("roleArn")
        sources = self._get_param("sources")
        tags = self._get_param("tags")
        input_type = self._get_param("type")
        a_input = self.medialive_backend.create_input(
            destinations=destinations,
            input_devices=input_devices,
            input_security_groups=input_security_groups,
            media_connect_flows=media_connect_flows,
            name=name,
            role_arn=role_arn,
            sources=sources,
            tags=tags,
            input_type=input_type,
        )
        return json.dumps({"input": a_input.to_dict()})

    def describe_input(self) -> str:
        input_id = self._get_param("inputId")
        a_input = self.medialive_backend.describe_input(input_id=input_id)
        return json.dumps(a_input.to_dict())

    def list_inputs(self) -> str:
        max_results = self._get_int_param("maxResults")
        inputs = self.medialive_backend.list_inputs(max_results=max_results)

        return json.dumps(dict(inputs=inputs, nextToken=None))

    def delete_input(self) -> str:
        input_id = self._get_param("inputId")
        self.medialive_backend.delete_input(input_id=input_id)
        return json.dumps({})

    def update_input(self) -> str:
        destinations = self._get_param("destinations")
        input_devices = self._get_param("inputDevices")
        input_id = self._get_param("inputId")
        input_security_groups = self._get_param("inputSecurityGroups")
        media_connect_flows = self._get_param("mediaConnectFlows")
        name = self._get_param("name")
        role_arn = self._get_param("roleArn")
        sources = self._get_param("sources")
        a_input = self.medialive_backend.update_input(
            destinations=destinations,
            input_devices=input_devices,
            input_id=input_id,
            input_security_groups=input_security_groups,
            media_connect_flows=media_connect_flows,
            name=name,
            role_arn=role_arn,
            sources=sources,
        )
        return json.dumps(dict(input=a_input.to_dict()))
