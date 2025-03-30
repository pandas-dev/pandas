import weakref
from typing import Dict, List

from moto.core.base_backend import BackendDict, BaseBackend
from moto.core.common_models import BaseModel

from .exceptions import (
    ApplicationNotFound,
    InvalidParameterValueError,
    ResourceNotFoundException,
)
from .utils import make_arn


class FakeEnvironment(BaseModel):
    def __init__(
        self,
        application: "FakeApplication",
        environment_name: str,
        solution_stack_name: str,
        tags: Dict[str, str],
    ):
        self.application = weakref.proxy(
            application
        )  # weakref to break circular dependencies
        self.environment_name = environment_name
        self.solution_stack_name = solution_stack_name
        self.tags = tags

    @property
    def application_name(self) -> str:
        return self.application.application_name

    @property
    def environment_arn(self) -> str:
        resource_path = f"{self.application_name}/{self.environment_name}"
        return make_arn(
            self.region, self.application.account_id, "environment", resource_path
        )

    @property
    def platform_arn(self) -> str:
        return "TODO"  # TODO

    @property
    def region(self) -> str:
        return self.application.region


class FakeApplication(BaseModel):
    def __init__(
        self,
        backend: "EBBackend",
        application_name: str,
    ):
        self.backend = weakref.proxy(backend)  # weakref to break cycles
        self.application_name = application_name
        self.environments: Dict[str, FakeEnvironment] = dict()
        self.account_id = self.backend.account_id
        self.region = self.backend.region_name
        self.arn = make_arn(
            self.region, self.account_id, "application", self.application_name
        )

    def create_environment(
        self, environment_name: str, solution_stack_name: str, tags: Dict[str, str]
    ) -> FakeEnvironment:
        if environment_name in self.environments:
            raise InvalidParameterValueError(message="")

        env = FakeEnvironment(
            application=self,
            environment_name=environment_name,
            solution_stack_name=solution_stack_name,
            tags=tags,
        )
        self.environments[environment_name] = env

        return env


class EBBackend(BaseBackend):
    def __init__(self, region_name: str, account_id: str):
        super().__init__(region_name, account_id)
        self.applications: Dict[str, FakeApplication] = dict()

    def create_application(self, application_name: str) -> FakeApplication:
        if application_name in self.applications:
            raise InvalidParameterValueError(
                f"Application {application_name} already exists."
            )
        new_app = FakeApplication(backend=self, application_name=application_name)
        self.applications[application_name] = new_app
        return new_app

    def create_environment(
        self,
        app: FakeApplication,
        environment_name: str,
        stack_name: str,
        tags: Dict[str, str],
    ) -> FakeEnvironment:
        return app.create_environment(
            environment_name=environment_name, solution_stack_name=stack_name, tags=tags
        )

    def describe_environments(self) -> List[FakeEnvironment]:
        envs = []
        for app in self.applications.values():
            for env in app.environments.values():
                envs.append(env)
        return envs

    def list_available_solution_stacks(self) -> None:
        # Implemented in response.py
        pass

    def update_tags_for_resource(
        self, resource_arn: str, tags_to_add: Dict[str, str], tags_to_remove: List[str]
    ) -> None:
        try:
            res = self._find_environment_by_arn(resource_arn)
        except KeyError:
            raise ResourceNotFoundException(
                f"Resource not found for ARN '{resource_arn}'."
            )

        for key, value in tags_to_add.items():
            res.tags[key] = value

        for key in tags_to_remove:
            del res.tags[key]

    def list_tags_for_resource(self, resource_arn: str) -> Dict[str, str]:
        try:
            res = self._find_environment_by_arn(resource_arn)
        except KeyError:
            raise ResourceNotFoundException(
                f"Resource not found for ARN '{resource_arn}'."
            )
        return res.tags

    def _find_environment_by_arn(self, arn: str) -> FakeEnvironment:
        for app in self.applications.keys():
            for env in self.applications[app].environments.values():
                if env.environment_arn == arn:
                    return env
        raise KeyError()

    def delete_application(
        self,
        application_name: str,
    ) -> None:
        if application_name:
            if application_name in self.applications:
                self.applications.pop(application_name)
            else:
                raise ApplicationNotFound(application_name)


eb_backends = BackendDict(EBBackend, "elasticbeanstalk")
