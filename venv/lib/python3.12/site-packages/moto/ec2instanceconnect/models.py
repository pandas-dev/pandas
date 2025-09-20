import json

from moto.core.base_backend import BackendDict, BaseBackend


class Ec2InstanceConnectBackend(BaseBackend):
    def send_ssh_public_key(self) -> str:
        return json.dumps(
            {"RequestId": "example-2a47-4c91-9700-e37e85162cb6", "Success": True}
        )


ec2instanceconnect_backends = BackendDict(Ec2InstanceConnectBackend, "ec2")
