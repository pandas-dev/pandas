from moto.core.responses import BaseResponse

from .models import Ec2InstanceConnectBackend, ec2instanceconnect_backends


class Ec2InstanceConnectResponse(BaseResponse):
    def __init__(self) -> None:
        super().__init__(service_name="ec2-instanceconnect")

    @property
    def ec2instanceconnect_backend(self) -> Ec2InstanceConnectBackend:
        return ec2instanceconnect_backends[self.current_account][self.region]

    def send_ssh_public_key(self) -> str:
        return self.ec2instanceconnect_backend.send_ssh_public_key()
