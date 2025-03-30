from moto.ec2.utils import utc_date_and_time

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

    def get_password_data(self) -> str:
        instance_id = self._get_param("InstanceId")
        password_data = self.ec2_backend.get_password_data(instance_id)
        template = self.response_template(GET_PASSWORD_DATA_RESPONSE)
        return template.render(
            password_data=password_data,
            instance_id=instance_id,
            timestamp=utc_date_and_time(),
        )


GET_PASSWORD_DATA_RESPONSE = """
<?xml version="1.0" encoding="UTF-8"?>
<GetPasswordDataResponse xmlns="http://ec2.amazonaws.com/doc/2016-11-15/">
  <requestId>6b9528d5-6818-4bd0-8936-cdedaEXAMPLE</requestId>
  <instanceId>{{ instance_id }}</instanceId>
  <timestamp>{{ timestamp }}</timestamp>
  <passwordData>{{ password_data }}</passwordData>
</GetPasswordDataResponse>
"""
