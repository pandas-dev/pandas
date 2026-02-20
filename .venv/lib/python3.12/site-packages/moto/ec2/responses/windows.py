from moto.core.responses import ActionResult
from moto.core.utils import utcnow

from ._base_response import EC2BaseResponse


class Windows(EC2BaseResponse):
    def bundle_instance(self) -> str:
        raise NotImplementedError("Windows.bundle_instance is not yet implemented")

    def cancel_bundle_task(self) -> str:
        raise NotImplementedError("Windows.cancel_bundle_task is not yet implemented")

    def describe_bundle_tasks(self) -> str:
        raise NotImplementedError(
            "Windows.describe_bundle_tasks is not yet implemented"
        )

    def get_password_data(self) -> ActionResult:
        instance_id = self._get_param("InstanceId")
        password_data = self.ec2_backend.get_password_data(instance_id)
        result = {
            "InstanceId": instance_id,
            "Timestamp": utcnow(),
            "PasswordData": password_data,
        }
        return ActionResult(result)
