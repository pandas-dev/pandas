from moto.core.common_models import BaseModel
from moto.moto_api._internal import mock_random as random


class WindowsBackend(BaseModel):
    def get_password_data(self, instance_id: str) -> str:
        instance = self.get_instance(instance_id)  # type: ignore[attr-defined]
        if instance.platform == "windows":
            return random.get_random_string(length=128)
        return ""
