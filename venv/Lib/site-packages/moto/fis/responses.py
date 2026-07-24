import json

from moto.core.responses import BaseResponse

from .models import FISBackend, fis_backends


class FISResponse(BaseResponse):
    def __init__(self) -> None:
        super().__init__(service_name="fis")

    @property
    def fis_backend(self) -> FISBackend:
        """Return backend instance specific for this region."""
        return fis_backends[self.current_account][self.region]

    def create_experiment_template(self) -> str:
        experiment_template = self.fis_backend.create_experiment_template(
            client_token=self._get_param("clientToken"),
            description=self._get_param("description"),
            stop_conditions=self._get_param("stopConditions"),
            targets=self._get_param("targets"),
            actions=self._get_param("actions"),
            role_arn=self._get_param("roleArn"),
            tags=self._get_param("tags"),
            log_configuration=self._get_param("logConfiguration"),
            experiment_options=self._get_param("experimentOptions"),
            experiment_report_configuration=self._get_param(
                "experimentReportConfiguration"
            ),
        )
        return json.dumps({"experimentTemplate": experiment_template})

    def delete_experiment_template(self) -> str:
        experiment_template = self.fis_backend.delete_experiment_template(
            id=self._get_param("id")
        )
        return json.dumps({"experimentTemplate": experiment_template})

    def tag_resource(self) -> str:
        self.fis_backend.tag_resource(
            resource_arn=self._get_param("resourceArn"),
            tags=self._get_param("tags"),
        )
        return "{}"

    def untag_resource(self) -> str:
        self.fis_backend.untag_resource(
            resource_arn=self._get_param("resourceArn"),
            tag_keys=self.querystring.get("tagKeys") or [],
        )
        return "{}"

    def list_tags_for_resource(self) -> str:
        tags = self.fis_backend.list_tags_for_resource(
            resource_arn=self._get_param("resourceArn")
        )
        return json.dumps({"tags": tags})
