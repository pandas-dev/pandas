from moto.core.responses import ActionResult, BaseResponse, EmptyResult
from moto.core.utils import tags_from_query_string

from .exceptions import InvalidParameterValueError
from .models import EBBackend, eb_backends

# AWS Solution Stack Details as of 2025-07-13
SOLUTION_STACK_DETAILS = [
    {
        "SolutionStackName": "64bit Windows Server 2019 v2.19.2 running IIS 10.0",
        "PermittedFileTypes": ["zip"],
    },
    {
        "SolutionStackName": "64bit Windows Server Core 2019 v2.19.2 running IIS 10.0",
        "PermittedFileTypes": ["zip"],
    },
    {
        "SolutionStackName": "64bit Windows Server 2025 v2.19.2 running IIS 10.0",
        "PermittedFileTypes": ["zip"],
    },
    {
        "SolutionStackName": "64bit Windows Server Core 2016 v2.19.2 running IIS 10.0",
        "PermittedFileTypes": ["zip"],
    },
    {
        "SolutionStackName": "64bit Windows Server Core 2025 v2.19.2 running IIS 10.0",
        "PermittedFileTypes": ["zip"],
    },
    {
        "SolutionStackName": "64bit Windows Server 2022 v2.19.2 running IIS 10.0",
        "PermittedFileTypes": ["zip"],
    },
    {
        "SolutionStackName": "64bit Windows Server 2016 v2.19.2 running IIS 10.0",
        "PermittedFileTypes": ["zip"],
    },
    {
        "SolutionStackName": "64bit Windows Server Core 2022 v2.19.2 running IIS 10.0",
        "PermittedFileTypes": ["zip"],
    },
    {
        "SolutionStackName": "64bit Amazon Linux 2017.03 v2.5.4 running Java 8",
        "PermittedFileTypes": ["jar", "zip"],
    },
    {
        "SolutionStackName": "64bit Amazon Linux 2023 v5.7.0 running Tomcat 9 Corretto 17",
        "PermittedFileTypes": ["war", "zip"],
    },
    {
        "SolutionStackName": "64bit Amazon Linux 2023 v5.7.0 running Tomcat 11 Corretto 17",
        "PermittedFileTypes": ["war", "zip"],
    },
    {
        "SolutionStackName": "64bit Amazon Linux 2023 v5.7.0 running Tomcat 10 Corretto 21",
        "PermittedFileTypes": ["war", "zip"],
    },
    {
        "SolutionStackName": "64bit Amazon Linux 2023 v5.7.0 running Tomcat 11 Corretto 21",
        "PermittedFileTypes": ["war", "zip"],
    },
    {
        "SolutionStackName": "64bit Amazon Linux 2023 v5.7.0 running Tomcat 10 Corretto 17",
        "PermittedFileTypes": ["war", "zip"],
    },
    {
        "SolutionStackName": "64bit Amazon Linux 2023 v5.7.0 running Tomcat 9 Corretto 11",
        "PermittedFileTypes": ["war", "zip"],
    },
    {
        "SolutionStackName": "64bit Amazon Linux 2023 v4.6.0 running Python 3.11",
        "PermittedFileTypes": ["zip"],
    },
    {
        "SolutionStackName": "64bit Amazon Linux 2023 v4.7.0 running PHP 8.2",
        "PermittedFileTypes": ["zip"],
    },
    {
        "SolutionStackName": "64bit Amazon Linux 2023 v4.7.0 running PHP 8.4",
        "PermittedFileTypes": ["zip"],
    },
    {
        "SolutionStackName": "64bit Amazon Linux 2023 v4.6.0 running Python 3.9",
        "PermittedFileTypes": ["zip"],
    },
    {
        "SolutionStackName": "64bit Amazon Linux 2023 v4.6.0 running Python 3.12",
        "PermittedFileTypes": ["zip"],
    },
    {
        "SolutionStackName": "64bit Amazon Linux 2023 v4.6.0 running Ruby 3.3",
        "PermittedFileTypes": ["zip"],
    },
    {
        "SolutionStackName": "64bit Amazon Linux 2023 v4.6.0 running Ruby 3.4",
        "PermittedFileTypes": ["zip"],
    },
    {
        "SolutionStackName": "64bit Amazon Linux 2023 v4.6.0 running Python 3.13",
        "PermittedFileTypes": ["zip"],
    },
    {
        "SolutionStackName": "64bit Amazon Linux 2023 v4.7.0 running PHP 8.3",
        "PermittedFileTypes": ["zip"],
    },
    {
        "SolutionStackName": "64bit Amazon Linux 2023 v4.6.0 running Ruby 3.2",
        "PermittedFileTypes": ["zip"],
    },
    {
        "SolutionStackName": "64bit Amazon Linux 2023 v6.6.0 running Node.js 22",
        "PermittedFileTypes": ["zip"],
    },
    {
        "SolutionStackName": "64bit Amazon Linux 2023 v4.7.0 running PHP 8.1",
        "PermittedFileTypes": ["zip"],
    },
    {
        "SolutionStackName": "64bit Amazon Linux 2023 v6.6.0 running Node.js 20",
        "PermittedFileTypes": ["zip"],
    },
    {
        "SolutionStackName": "64bit Amazon Linux 2023 v4.6.0 running Corretto 17",
        "PermittedFileTypes": ["jar", "zip"],
    },
    {
        "SolutionStackName": "64bit Amazon Linux 2023 v3.5.0 running .NET 9",
        "PermittedFileTypes": ["zip"],
    },
    {
        "SolutionStackName": "64bit Amazon Linux 2023 v6.6.0 running Node.js 18",
        "PermittedFileTypes": ["zip"],
    },
    {
        "SolutionStackName": "64bit Amazon Linux 2023 v4.6.0 running Corretto 8",
        "PermittedFileTypes": ["jar", "zip"],
    },
    {
        "SolutionStackName": "64bit Amazon Linux 2023 v3.5.0 running .NET 8",
        "PermittedFileTypes": ["zip"],
    },
    {
        "SolutionStackName": "64bit Amazon Linux 2023 v4.6.0 running Corretto 11",
        "PermittedFileTypes": ["jar", "zip"],
    },
    {
        "SolutionStackName": "64bit Amazon Linux 2023 v4.6.0 running Docker",
        "PermittedFileTypes": [],
    },
    {
        "SolutionStackName": "64bit Amazon Linux 2023 v4.2.0 running ECS",
        "PermittedFileTypes": [],
    },
    {
        "SolutionStackName": "64bit Amazon Linux 2023 v4.4.0 running Go 1",
        "PermittedFileTypes": ["zip"],
    },
    {
        "SolutionStackName": "64bit Amazon Linux 2023 v4.6.0 running Corretto 21",
        "PermittedFileTypes": ["jar", "zip"],
    },
    {
        "SolutionStackName": "64bit Amazon Linux 2 v4.9.0 running Tomcat 9 Corretto 11",
        "PermittedFileTypes": ["war", "zip"],
    },
    {
        "SolutionStackName": "64bit Amazon Linux 2 v2.11.0 running .NET Core",
        "PermittedFileTypes": ["zip"],
    },
    {
        "SolutionStackName": "64bit Amazon Linux 2 v3.9.0 running Corretto 11",
        "PermittedFileTypes": ["jar", "zip"],
    },
    {
        "SolutionStackName": "64bit Amazon Linux 2 v3.9.0 running Corretto 8",
        "PermittedFileTypes": ["jar", "zip"],
    },
    {
        "SolutionStackName": "64bit Amazon Linux 2 v3.13.0 running Go 1",
        "PermittedFileTypes": ["zip"],
    },
    {
        "SolutionStackName": "64bit Amazon Linux 2 v5.11.0 running Node.js 18",
        "PermittedFileTypes": ["zip"],
    },
    {
        "SolutionStackName": "64bit Amazon Linux 2 v3.5.0 running ECS",
        "PermittedFileTypes": [],
    },
    {
        "SolutionStackName": "64bit Amazon Linux 2 v3.9.0 running Corretto 17",
        "PermittedFileTypes": ["jar", "zip"],
    },
    {
        "SolutionStackName": "64bit Amazon Linux 2 v3.10.0 running PHP 8.1",
        "PermittedFileTypes": ["zip"],
    },
    {
        "SolutionStackName": "64bit Amazon Linux 2 v4.9.0 running Tomcat 9 Corretto 8",
        "PermittedFileTypes": ["war", "zip"],
    },
    {
        "SolutionStackName": "64bit Amazon Linux 2 v4.2.0 running Docker",
        "PermittedFileTypes": [],
    },
    {
        "SolutionStackName": "64bit Amazon Linux 2 v3.3.14 running Python 3.8",
        "PermittedFileTypes": ["zip"],
    },
    {
        "SolutionStackName": "64bit Amazon Linux 2 v5.0.2 running Node.js 12",
        "PermittedFileTypes": ["zip"],
    },
    {
        "SolutionStackName": "64bit Amazon Linux 2 v5.0.2 running Node.js 10",
        "PermittedFileTypes": ["zip"],
    },
    {
        "SolutionStackName": "64bit Amazon Linux 2 v0.1.0 running Node.js 12 (BETA)",
        "PermittedFileTypes": ["zip"],
    },
    {
        "SolutionStackName": "64bit Amazon Linux 2018.03 v2.6.33 running Packer 1.0.3",
        "PermittedFileTypes": [],
    },
]


class EBResponse(BaseResponse):
    def __init__(self) -> None:
        super().__init__(service_name="elasticbeanstalk")

    @property
    def elasticbeanstalk_backend(self) -> EBBackend:
        """
        :rtype: EBBackend
        """
        return eb_backends[self.current_account][self.region]

    def create_application(self) -> ActionResult:
        app = self.elasticbeanstalk_backend.create_application(
            application_name=self._get_param("ApplicationName")
        )

        result = {"Application": app}
        return ActionResult(result)

    def describe_applications(self) -> ActionResult:
        applications = self.elasticbeanstalk_backend.applications.values()
        result = {"Applications": applications}
        return ActionResult(result)

    def create_environment(self) -> ActionResult:
        application_name = self._get_param("ApplicationName")
        try:
            app = self.elasticbeanstalk_backend.applications[application_name]
        except KeyError:
            raise InvalidParameterValueError(
                f"No Application named '{application_name}' found."
            )

        tags = tags_from_query_string(self.querystring, prefix="Tags.member")
        env = self.elasticbeanstalk_backend.create_environment(
            app,
            environment_name=self._get_param("EnvironmentName"),
            stack_name=self._get_param("SolutionStackName"),
            tags=tags,
        )

        result = env
        return ActionResult(result)

    def describe_environments(self) -> ActionResult:
        envs = self.elasticbeanstalk_backend.describe_environments()

        result = {"Environments": envs}
        return ActionResult(result)

    def list_available_solution_stacks(self) -> ActionResult:
        solution_stacks = {
            "SolutionStacks": [
                stack["SolutionStackName"] for stack in SOLUTION_STACK_DETAILS
            ],
            "SolutionStackDetails": SOLUTION_STACK_DETAILS,
        }
        return ActionResult(solution_stacks)

    def update_tags_for_resource(self) -> ActionResult:
        resource_arn = self._get_param("ResourceArn")
        tags_to_add = tags_from_query_string(
            self.querystring, prefix="TagsToAdd.member"
        )
        tags_to_remove = self._get_multi_param("TagsToRemove.member")
        self.elasticbeanstalk_backend.update_tags_for_resource(
            resource_arn, tags_to_add, tags_to_remove
        )
        return EmptyResult()

    def list_tags_for_resource(self) -> ActionResult:
        resource_arn = self._get_param("ResourceArn")
        tags = self.elasticbeanstalk_backend.list_tags_for_resource(resource_arn)
        result = {
            "ResourceArn": resource_arn,
            "ResourceTags": [{"Key": k, "Value": v} for k, v in tags.items()],
        }
        return ActionResult(result)

    def delete_application(self) -> ActionResult:
        application_name = self._get_param("ApplicationName")
        self.elasticbeanstalk_backend.delete_application(
            application_name=application_name,
        )
        return EmptyResult()
