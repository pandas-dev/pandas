from moto.core.responses import ActionResult, BaseResponse, EmptyResult, PaginatedResult

from .models import MediaLiveBackend, medialive_backends


class MediaLiveResponse(BaseResponse):
    def __init__(self) -> None:
        super().__init__(service_name="medialive")
        self.automated_parameter_parsing = True

    @property
    def medialive_backend(self) -> MediaLiveBackend:
        return medialive_backends[self.current_account][self.region]

    def create_channel(self) -> ActionResult:
        cdi_input_specification = self._get_param("CdiInputSpecification")
        channel_class = self._get_param("ChannelClass")
        destinations = self._get_param("Destinations")
        encoder_settings = self._get_param("EncoderSettings")
        input_attachments = self._get_param("InputAttachments")
        input_specification = self._get_param("InputSpecification")
        log_level = self._get_param("LogLevel")
        name = self._get_param("Name")
        role_arn = self._get_param("RoleArn")
        tags = self._get_param("Tags")
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

        return ActionResult({"Channel": channel})

    def list_channels(self) -> ActionResult:
        channels = self.medialive_backend.list_channels()
        return PaginatedResult({"Channels": channels})

    def describe_channel(self) -> ActionResult:
        channel_id = self._get_param("ChannelId")
        channel = self.medialive_backend.describe_channel(channel_id=channel_id)
        return ActionResult(channel)

    def delete_channel(self) -> ActionResult:
        channel_id = self._get_param("ChannelId")
        channel = self.medialive_backend.delete_channel(channel_id=channel_id)
        return ActionResult(channel)

    def start_channel(self) -> ActionResult:
        channel_id = self._get_param("ChannelId")
        channel = self.medialive_backend.start_channel(channel_id=channel_id)
        return ActionResult(channel)

    def stop_channel(self) -> ActionResult:
        channel_id = self._get_param("ChannelId")
        channel = self.medialive_backend.stop_channel(channel_id=channel_id)
        return ActionResult(channel)

    def update_channel(self) -> ActionResult:
        channel_id = self._get_param("ChannelId")
        cdi_input_specification = self._get_param("CdiInputSpecification")
        destinations = self._get_param("Destinations")
        encoder_settings = self._get_param("EncoderSettings")
        input_attachments = self._get_param("InputAttachments")
        input_specification = self._get_param("InputSpecification")
        log_level = self._get_param("LogLevel")
        name = self._get_param("Name")
        role_arn = self._get_param("RoleArn")
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
        return ActionResult({"Channel": channel})

    def create_input(self) -> ActionResult:
        destinations = self._get_param("Destinations")
        input_devices = self._get_param("InputDevices")
        input_security_groups = self._get_param("InputSecurityGroups")
        media_connect_flows = self._get_param("MediaConnectFlows")
        name = self._get_param("Name")
        role_arn = self._get_param("RoleArn")
        sources = self._get_param("Sources")
        tags = self._get_param("Tags")
        input_type = self._get_param("Type")
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
        return ActionResult({"Input": a_input})

    def describe_input(self) -> ActionResult:
        input_id = self._get_param("InputId")
        a_input = self.medialive_backend.describe_input(input_id=input_id)
        return ActionResult(a_input)

    def list_inputs(self) -> ActionResult:
        inputs = self.medialive_backend.list_inputs()
        return PaginatedResult({"Inputs": inputs})

    def delete_input(self) -> ActionResult:
        input_id = self._get_param("InputId")
        self.medialive_backend.delete_input(input_id=input_id)
        return EmptyResult()

    def update_input(self) -> ActionResult:
        destinations = self._get_param("Destinations")
        input_devices = self._get_param("InputDevices")
        input_id = self._get_param("InputId")
        input_security_groups = self._get_param("InputSecurityGroups")
        media_connect_flows = self._get_param("MediaConnectFlows")
        name = self._get_param("Name")
        role_arn = self._get_param("RoleArn")
        sources = self._get_param("Sources")
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
        return ActionResult({"Input": a_input})
