from collections import OrderedDict
from typing import Any, Optional

from moto.core.common_models import CloudFormationModel
from moto.core.types import Base64EncodedString

from ..exceptions import (
    InvalidLaunchTemplateIdNotFound,
    InvalidLaunchTemplateNameAlreadyExistsError,
    InvalidLaunchTemplateNameNotFoundError,
    InvalidLaunchTemplateNameNotFoundWithNameError,
    InvalidLaunchTemplateVersionNotFound,
    InvalidParameterValue,
    MissingParameter,
    MissingParameterError,
)
from ..utils import (
    convert_tag_spec,
    generic_filter,
    parse_user_data,
    random_launch_template_id,
    random_launch_template_name,
    utc_date_and_time,
)
from .core import TaggedEC2Resource


class LaunchTemplateVersion:
    def __init__(
        self,
        template: "LaunchTemplate",
        number: int,
        data: dict[str, Any],
        description: str,
    ):
        self.template = template
        self.number = number
        self.data = data
        self.description = description
        self.create_time = utc_date_and_time()
        self.instance_tags = convert_tag_spec(data.get("TagSpecifications", [])).get(
            "instance", {}
        )

    @property
    def image_id(self) -> str:
        return self.data.get("ImageId", "")

    @property
    def instance_type(self) -> str:
        return self.data.get("InstanceType", "")

    @property
    def security_groups(self) -> list[str]:
        return self.data.get("SecurityGroups", [])

    @property
    def user_data(self) -> Optional[Base64EncodedString]:
        user_data = self.data.get("UserData")
        # UserData can be specified via multiple services/api endpoints,
        # so we make an assertion here that it's in the format we expect.
        assert user_data is None or isinstance(user_data, Base64EncodedString)
        return user_data


class LaunchTemplate(TaggedEC2Resource, CloudFormationModel):
    def __init__(
        self,
        backend: Any,
        name: str,
        template_data: dict[str, Any],
        version_description: str,
        tag_spec: dict[str, dict[str, str]],
    ):
        self.ec2_backend = backend
        self.name = name
        self.id = random_launch_template_id()
        self.create_time = utc_date_and_time()
        tag_map: dict[str, str] = tag_spec.get("launch-template", {})
        self.add_tags(tag_map)
        self.tags = self.get_tags()

        self.versions: list[LaunchTemplateVersion] = []
        self.create_version(template_data, version_description)
        self.default_version_number = 1

    def create_version(
        self, data: dict[str, Any], description: str
    ) -> LaunchTemplateVersion:
        num = len(self.versions) + 1
        version = LaunchTemplateVersion(self, num, data, description)
        self.versions.append(version)
        return version

    def is_default(self, version: LaunchTemplateVersion) -> bool:
        return self.default_version_number == version.number

    def get_version(self, num: Any) -> LaunchTemplateVersion:
        if str(num).lower() == "$latest":
            return self.versions[-1]
        if str(num).lower() == "$default":
            return self.default_version()
        try:
            return self.versions[int(num) - 1]
        except IndexError:
            raise InvalidLaunchTemplateVersionNotFound(template_id=self.id, version=num)

    def default_version(self) -> LaunchTemplateVersion:
        return self.versions[self.default_version_number - 1]

    def latest_version(self) -> LaunchTemplateVersion:
        return self.versions[-1]

    @property
    def latest_version_number(self) -> int:
        return self.latest_version().number

    @property
    def physical_resource_id(self) -> str:
        return self.id

    def get_filter_value(
        self, filter_name: str, method_name: Optional[str] = None
    ) -> Any:
        if filter_name == "launch-template-name":
            return self.name
        else:
            return super().get_filter_value(filter_name, "DescribeLaunchTemplates")

    @staticmethod
    def cloudformation_name_type() -> str:
        return "LaunchTemplateName"

    @staticmethod
    def cloudformation_type() -> str:
        # https://docs.aws.amazon.com/AWSCloudFormation/latest/UserGuide/aws-resource-ec2-launchtemplate.html
        return "AWS::EC2::LaunchTemplate"

    @classmethod
    def create_from_cloudformation_json(  # type: ignore[misc]
        cls,
        resource_name: str,
        cloudformation_json: Any,
        account_id: str,
        region_name: str,
        **kwargs: Any,
    ) -> "LaunchTemplate":
        from ..models import ec2_backends

        backend = ec2_backends[account_id][region_name]

        properties = cloudformation_json["Properties"]
        name = properties.get("LaunchTemplateName")
        data = properties.get("LaunchTemplateData")
        data["UserData"] = parse_user_data(data.get("UserData"))
        description = properties.get("VersionDescription")
        tag_spec = convert_tag_spec(
            properties.get("TagSpecifications", {}), tag_key="Tags"
        )

        if name is None:
            name = random_launch_template_name()

        launch_template = backend.create_launch_template(
            name, description, data, tag_spec
        )

        return launch_template

    @classmethod
    def update_from_cloudformation_json(  # type: ignore[misc]
        cls,
        original_resource: Any,
        new_resource_name: str,
        cloudformation_json: Any,
        account_id: str,
        region_name: str,
    ) -> "LaunchTemplate":
        from ..models import ec2_backends

        backend = ec2_backends[account_id][region_name]

        properties = cloudformation_json["Properties"]

        data = properties.get("LaunchTemplateData")
        data["UserData"] = parse_user_data(data.get("UserData"))
        description = properties.get("VersionDescription")

        launch_template = backend.get_launch_template(original_resource.id)

        launch_template.create_version(data, description)

        return launch_template

    def delete(
        self,
        account_id: str,
        region_name: str,
    ) -> None:
        from ..models import ec2_backends

        backend = ec2_backends[account_id][region_name]
        backend.delete_launch_template(self.name, None)

    @classmethod
    def has_cfn_attr(cls, attr: str) -> bool:
        return attr in [
            "DefaultVersionNumber",
            "LaunchTemplateId",
            "LatestVersionNumber",
        ]

    def get_cfn_attribute(self, attribute_name: str) -> str:
        from moto.cloudformation.exceptions import UnformattedGetAttTemplateException

        if attribute_name == "DefaultVersionNumber":
            return str(self.default_version_number)
        if attribute_name == "LaunchTemplateId":
            return self.id
        if attribute_name == "LatestVersionNumber":
            return str(self.latest_version_number)
        raise UnformattedGetAttTemplateException()


class LaunchTemplateBackend:
    def __init__(self) -> None:
        self.launch_template_name_to_ids: dict[str, str] = {}
        self.launch_templates: dict[str, LaunchTemplate] = OrderedDict()
        self.launch_template_insert_order: list[str] = []

    def create_launch_template(
        self,
        name: str,
        description: str,
        template_data: dict[str, Any],
        tag_spec: dict[str, Any],
    ) -> LaunchTemplate:
        if name in self.launch_template_name_to_ids:
            raise InvalidLaunchTemplateNameAlreadyExistsError()
        template = LaunchTemplate(self, name, template_data, description, tag_spec)
        self.launch_templates[template.id] = template
        self.launch_template_name_to_ids[template.name] = template.id
        self.launch_template_insert_order.append(template.id)
        return template

    def modify_launch_template(
        self,
        default_version: str,
        template_name: Optional[str] = None,
        template_id: Optional[str] = None,
    ) -> LaunchTemplate:
        if template_name:
            template_id = self.launch_template_name_to_ids.get(template_name)
        if template_id is None:
            raise MissingParameterError("launch template ID or launch template name")
        if template_id not in self.launch_templates:
            raise InvalidLaunchTemplateNameNotFoundError()
        template = self.launch_templates[template_id]
        if default_version not in [str(v.number) for v in template.versions]:
            raise InvalidLaunchTemplateVersionNotFound(version=default_version)

        template.default_version_number = int(default_version)
        template = self.launch_templates[template_id]
        return template

    def get_launch_template(self, template_id: str) -> LaunchTemplate:
        if template_id not in self.launch_templates:
            raise InvalidLaunchTemplateIdNotFound(template_id)
        return self.launch_templates[template_id]

    def get_launch_template_by_name(self, name: str) -> LaunchTemplate:
        if name not in self.launch_template_name_to_ids:
            raise InvalidLaunchTemplateNameNotFoundWithNameError(name)
        return self.get_launch_template(self.launch_template_name_to_ids[name])

    def delete_launch_template(self, name: str, tid: Optional[str]) -> LaunchTemplate:
        if name:
            tid = self.launch_template_name_to_ids.get(name)
        if tid is None:
            raise MissingParameterError("launch template ID or launch template name")
        if tid not in self.launch_templates:
            raise InvalidLaunchTemplateNameNotFoundError()
        template = self.launch_templates.pop(tid)
        self.launch_template_name_to_ids.pop(template.name, None)
        return template

    def describe_launch_templates(
        self,
        template_names: Optional[list[str]] = None,
        template_ids: Optional[list[str]] = None,
        filters: Any = None,
    ) -> list[LaunchTemplate]:
        if template_names and not template_ids:
            template_ids = []
            for name in template_names:
                if name not in self.launch_template_name_to_ids:
                    raise InvalidLaunchTemplateNameNotFoundError()
                template_ids.append(self.launch_template_name_to_ids[name])

        if template_ids:
            templates = [
                self.launch_templates[tid]
                for tid in template_ids
                if tid in self.launch_templates
            ]
        else:
            templates = list(self.launch_templates.values())

        return generic_filter(filters, templates)

    def describe_launch_template_versions(
        self,
        template_name: str,
        template_id: str,
        versions: list[str],
        min_version: int,
        max_version: int,
        max_results: int,
    ) -> list[tuple[LaunchTemplate, LaunchTemplateVersion]]:
        """
        The Filters-parameter is not yet implemented
        """
        wrong_param_msg = "To describe the launch template data for all your launch templates, for ‘--versions’ specify ‘$Latest’, ‘$Default’, or both, and omit ‘--launch-template-id’, ‘--launch-template-name’, and version numbers. To describe the launch template data for a specific launch template, specify ‘--launch-template-id’ or ‘--launch-template-name’, and for ‘--versions’ specify one or more of the following values: ‘$Latest’, ‘$Default’, or one or more version numbers."
        if template_name:
            template = self.get_launch_template_by_name(template_name)
        elif template_id:
            template = self.get_launch_template(template_id)
        elif not versions:
            raise MissingParameter(wrong_param_msg)
        elif versions and ("$Latest" not in versions and "$Default" not in versions):
            raise InvalidParameterValue(wrong_param_msg)
        else:
            template = None

        ret_versions: list[tuple[LaunchTemplate, LaunchTemplateVersion]] = []
        if versions and template:
            for v in versions:
                ret_versions.append((template, template.get_version(v)))
        elif not template:
            # Version has to have either $Latest and/or $Default at this point. This was already validated earlier, so here we can just iterate over both
            for template in self.launch_templates.values():
                if "$Latest" in versions:
                    ret_versions.append((template, template.get_version("$Latest")))
                if "$Default" in versions:
                    ret_versions.append((template, template.get_version("$Default")))
        elif min_version:
            if max_version:
                vMax = max_version
            else:
                vMax = min_version + max_results

            vMin = min_version - 1
            ret_versions = [(template, ver) for ver in template.versions[vMin:vMax]]
        elif max_version:
            vMax = max_version
            ret_versions = [(template, ver) for ver in template.versions[:vMax]]
        elif template is not None:
            ret_versions = [(template, ver) for ver in template.versions]

        ret_versions = ret_versions[:max_results]
        return ret_versions

    def get_launch_template_data(self, instance_id: str) -> Any:
        return self.get_instance(instance_id)  # type: ignore[attr-defined]
