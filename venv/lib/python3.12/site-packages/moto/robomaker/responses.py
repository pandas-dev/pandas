import json

from moto.core.responses import BaseResponse

from .models import RoboMakerBackend, robomaker_backends


class RoboMakerResponse(BaseResponse):
    def __init__(self) -> None:
        super().__init__(service_name="robomaker")

    @property
    def robomaker_backend(self) -> RoboMakerBackend:
        return robomaker_backends[self.current_account][self.region]

    def create_robot_application(self) -> str:
        name = self._get_param("name")
        sources = self._get_param("sources")
        robot_software_suite = self._get_param("robotSoftwareSuite")
        app = self.robomaker_backend.create_robot_application(
            name=name,
            sources=sources,
            robot_software_suite=robot_software_suite,
        )
        return json.dumps(app.to_dict())

    def describe_robot_application(self) -> str:
        application = self._get_param("application")
        app = self.robomaker_backend.describe_robot_application(
            application=application,
        )
        return json.dumps(app.to_dict())

    def delete_robot_application(self) -> str:
        application = self._get_param("application")
        self.robomaker_backend.delete_robot_application(
            application=application,
        )
        return "{}"

    def list_robot_applications(self) -> str:
        robot_applications = self.robomaker_backend.list_robot_applications()
        return json.dumps(
            dict(
                robotApplicationSummaries=[app.to_dict() for app in robot_applications]
            )
        )
