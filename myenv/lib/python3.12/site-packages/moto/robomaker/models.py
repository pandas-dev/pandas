from typing import Any, Dict, Iterable, List

from moto.core.base_backend import BackendDict, BaseBackend
from moto.core.common_models import BaseModel
from moto.core.utils import unix_time
from moto.utilities.utils import get_partition


class RobotApplication(BaseModel):
    def __init__(
        self,
        account_id: str,
        region: str,
        name: str,
        sources: List[Dict[str, str]],
        robot_software_suite: Dict[str, str],
    ):
        self.account_id = account_id
        self.region = region
        self.name = name
        self.sources = sources
        self.robot_software_suite = robot_software_suite
        self.created_on = unix_time()

    @property
    def arn(self) -> str:
        return f"arn:{get_partition(self.region)}:robomaker:{self.region}:{self.account_id}:robot-application/{self.name}/{self.created_on}"

    def to_dict(self) -> Dict[str, Any]:
        return {
            "arn": self.arn,
            "name": self.name,
            "sources": self.sources,
            "robotSoftwareSuite": self.robot_software_suite,
        }


class RoboMakerBackend(BaseBackend):
    def __init__(self, region_name: str, account_id: str):
        super().__init__(region_name, account_id)
        self.robot_applications: Dict[str, RobotApplication] = {}

    def create_robot_application(
        self,
        name: str,
        sources: List[Dict[str, str]],
        robot_software_suite: Dict[str, str],
    ) -> RobotApplication:
        """
        The tags and environment parameters are not yet implemented
        """
        app = RobotApplication(
            account_id=self.account_id,
            region=self.region_name,
            name=name,
            sources=sources,
            robot_software_suite=robot_software_suite,
        )
        self.robot_applications[name] = app
        return app

    def describe_robot_application(self, application: str) -> RobotApplication:
        return self.robot_applications[application]

    def delete_robot_application(self, application: str) -> None:
        self.robot_applications.pop(application)

    def list_robot_applications(self) -> Iterable[RobotApplication]:
        """
        Currently returns all applications - none of the parameters are taken into account
        """
        return self.robot_applications.values()


robomaker_backends = BackendDict(RoboMakerBackend, "robomaker")
