from typing import Any, Dict, List, Optional

from moto.core.base_backend import BackendDict, BaseBackend
from moto.core.common_models import BaseModel
from moto.core.utils import unix_time
from moto.moto_api._internal import mock_random
from moto.utilities.paginator import paginate
from moto.utilities.tagging_service import TaggingService
from moto.utilities.utils import get_partition

from .exceptions import AppNotFound, AppVersionNotFound, ResiliencyPolicyNotFound

PAGINATION_MODEL = {
    "list_apps": {
        "input_token": "next_token",
        "limit_key": "max_results",
        "limit_default": 100,
        "unique_attribute": "arn",
    },
    "list_resiliency_policies": {
        "input_token": "next_token",
        "limit_key": "max_results",
        "limit_default": 100,
        "unique_attribute": "arn",
    },
}


class AppComponent(BaseModel):
    def __init__(self, _id: str, name: str, _type: str):
        self.id = _id
        self.name = name
        self.type = _type

    def to_json(self) -> Dict[str, Any]:
        return {
            "id": self.id,
            "name": self.name,
            "type": self.type,
        }


class App(BaseModel):
    def __init__(
        self,
        backend: "ResilienceHubBackend",
        assessment_schedule: str,
        description: str,
        event_subscriptions: List[Dict[str, Any]],
        name: str,
        permission_model: Dict[str, Any],
        policy_arn: str,
    ):
        self.backend = backend
        self.arn = f"arn:{get_partition(backend.region_name)}:resiliencehub:{backend.region_name}:{backend.account_id}:app/{mock_random.uuid4()}"
        self.assessment_schedule = assessment_schedule or "Disabled"
        self.compliance_status = "NotAssessed"
        self.description = description
        self.creation_time = unix_time()
        self.event_subscriptions = event_subscriptions
        self.name = name
        self.permission_model = permission_model
        self.policy_arn = policy_arn
        self.resilience_score = 0.0
        self.status = "Active"
        self.app_versions: List[AppVersion] = []

        app_version = AppVersion(app_arn=self.arn, version_name=None, identifier=0)
        self.app_versions.append(app_version)

    def get_version(self, version_name: str) -> "AppVersion":
        for v in self.app_versions:
            if v.app_version == version_name:
                return v
        raise AppVersionNotFound

    def to_json(self) -> Dict[str, Any]:
        resp = {
            "appArn": self.arn,
            "assessmentSchedule": self.assessment_schedule,
            "complianceStatus": self.compliance_status,
            "creationTime": self.creation_time,
            "driftStatus": "NotChecked",
            "name": self.name,
            "resilienceScore": self.resilience_score,
            "status": self.status,
            "tags": self.backend.list_tags_for_resource(self.arn),
        }
        if self.description is not None:
            resp["description"] = self.description
        if self.event_subscriptions:
            resp["eventSubscriptions"] = self.event_subscriptions
        if self.permission_model:
            resp["permissionModel"] = self.permission_model
        if self.policy_arn:
            resp["policyArn"] = self.policy_arn
        return resp


class Resource:
    def __init__(
        self,
        logical_resource_id: Dict[str, Any],
        physical_resource_id: str,
        resource_type: str,
        components: List[AppComponent],
    ):
        self.logical_resource_id = logical_resource_id
        self.physical_resource_id = physical_resource_id
        self.resource_type = resource_type
        self.components = components

    def to_json(self) -> Dict[str, Any]:
        return {
            "appComponents": [c.to_json() for c in self.components],
            "resourceType": self.resource_type,
            "logicalResourceId": self.logical_resource_id,
            "physicalResourceId": {"identifier": self.physical_resource_id},
            "resourceName": self.logical_resource_id["identifier"],
        }


class AppVersion(BaseModel):
    def __init__(self, app_arn: str, version_name: Optional[str], identifier: int):
        self.app_arn = app_arn
        self.eks_sources: List[Dict[str, Any]] = []
        self.source_arns: List[str] = []
        self.terraform_sources: List[Dict[str, str]] = []
        self.app_version = "release" if version_name else "draft"
        self.identifier = identifier
        self.creation_time = unix_time()
        self.version_name = version_name
        self.app_components: List[AppComponent] = []
        self.status = "Pending"
        self.resources: List[Resource] = []

    def to_json(self) -> Dict[str, Any]:
        resp = {
            "appVersion": self.app_version,
            "creationTime": self.creation_time,
            "identifier": self.identifier,
        }
        if self.version_name:
            resp["versionName"] = self.version_name
        return resp


class Policy(BaseModel):
    def __init__(
        self,
        backend: "ResilienceHubBackend",
        policy: Dict[str, Dict[str, int]],
        policy_name: str,
        data_location_constraint: str,
        policy_description: str,
        tier: str,
    ):
        self.arn = f"arn:{get_partition(backend.region_name)}:resiliencehub:{backend.region_name}:{backend.account_id}:resiliency-policy/{mock_random.uuid4()}"
        self.backend = backend
        self.data_location_constraint = data_location_constraint
        self.creation_time = unix_time()
        self.policy = policy
        self.policy_description = policy_description
        self.policy_name = policy_name
        self.tier = tier

    def to_json(self) -> Dict[str, Any]:
        resp = {
            "creationTime": self.creation_time,
            "policy": self.policy,
            "policyArn": self.arn,
            "policyName": self.policy_name,
            "tags": self.backend.list_tags_for_resource(self.arn),
            "tier": self.tier,
        }
        if self.data_location_constraint:
            resp["dataLocationConstraint"] = self.data_location_constraint
        if self.policy_description:
            resp["policyDescription"] = self.policy_description
        return resp


class ResilienceHubBackend(BaseBackend):
    def __init__(self, region_name: str, account_id: str):
        super().__init__(region_name, account_id)
        self.apps: Dict[str, App] = dict()
        self.policies: Dict[str, Policy] = dict()
        self.tagger = TaggingService()

        self.app_assessments_queue: List[List[Dict[str, Any]]] = []
        self.app_assessments_results: Dict[str, List[Dict[str, Any]]] = {}

    def create_app(
        self,
        assessment_schedule: str,
        description: str,
        event_subscriptions: List[Dict[str, Any]],
        name: str,
        permission_model: Dict[str, Any],
        policy_arn: str,
        tags: Dict[str, str],
    ) -> App:
        """
        The ClientToken-parameter is not yet implemented
        """
        app = App(
            backend=self,
            assessment_schedule=assessment_schedule,
            description=description,
            event_subscriptions=event_subscriptions,
            name=name,
            permission_model=permission_model,
            policy_arn=policy_arn,
        )
        self.apps[app.arn] = app
        self.tag_resource(app.arn, tags)
        return app

    def create_resiliency_policy(
        self,
        data_location_constraint: str,
        policy: Dict[str, Any],
        policy_description: str,
        policy_name: str,
        tags: Dict[str, str],
        tier: str,
    ) -> Policy:
        """
        The ClientToken-parameter is not yet implemented
        """
        pol = Policy(
            backend=self,
            data_location_constraint=data_location_constraint,
            policy=policy,
            policy_description=policy_description,
            policy_name=policy_name,
            tier=tier,
        )
        self.policies[pol.arn] = pol
        self.tag_resource(pol.arn, tags)
        return pol

    @paginate(PAGINATION_MODEL)
    def list_apps(self, app_arn: str, name: str, reverse_order: bool) -> List[App]:
        """
        The FromAssessmentTime/ToAssessmentTime-parameters are not yet implemented
        """
        if name:
            app_summaries = [a for a in self.apps.values() if a.name == name]
        elif app_arn:
            app_summaries = [self.apps[app_arn]]
        else:
            app_summaries = list(self.apps.values())
        if reverse_order:
            app_summaries.reverse()
        return app_summaries

    def list_app_assessments(self, request_identifier: str) -> List[Dict[str, Any]]:
        """
        Moto will not actually execute any assessments, so this operation will return an empty list by default.
        You can use a dedicated API to override this, by configuring a queue of expected results.

        A request to `list_app_assessments` will take the first result from that queue, with subsequent calls with the same parameters returning the same result.

        Calling `list_app_assessments` with a different set of parameters will return the second result from that queue - and so on, or an empty list of the queue is empty.

        Configure this queue by making an HTTP request to `/moto-api/static/resilience-hub-assessments/response`. An example invocation looks like this:

        .. sourcecode:: python

            summary1 = {"appArn": "app_arn1", "appVersion": "some version", ...}
            summary2 = {"appArn": "app_arn2", ...}
            results = {"results": [[summary1, summary2], [summary2]], "region": "us-east-1"}
            resp = requests.post(
                "http://motoapi.amazonaws.com/moto-api/static/resilience-hub-assessments/response",
                json=results,
            )

            assert resp.status_code == 201

            client = boto3.client("lambda", region_name="us-east-1")
            # First result
            resp = client.list_app_assessments() # [summary1, summary2]
            # Second result
            resp = client.list_app_assessments(assessmentStatus="Pending") # [summary2]

        If you're using MotoServer, make sure to make this request to where MotoServer is running:

        .. sourcecode:: python

            http://localhost:5000/moto-api/static/resilience-hub-assessments/response

        """
        if request_identifier in self.app_assessments_results:
            return self.app_assessments_results[request_identifier]
        if self.app_assessments_queue:
            self.app_assessments_results[request_identifier] = (
                self.app_assessments_queue.pop(0)
            )
            return self.app_assessments_results[request_identifier]
        return []

    def describe_app(self, app_arn: str) -> App:
        if app_arn not in self.apps:
            raise AppNotFound(app_arn)
        return self.apps[app_arn]

    @paginate(pagination_model=PAGINATION_MODEL)
    def list_resiliency_policies(self, policy_name: str) -> List[Policy]:
        if policy_name:
            return [p for p in self.policies.values() if p.policy_name == policy_name]
        return list(self.policies.values())

    def describe_resiliency_policy(self, policy_arn: str) -> Policy:
        if policy_arn not in self.policies:
            raise ResiliencyPolicyNotFound(policy_arn)
        return self.policies[policy_arn]

    def tag_resource(self, resource_arn: str, tags: Dict[str, str]) -> None:
        self.tagger.tag_resource(
            resource_arn, TaggingService.convert_dict_to_tags_input(tags)
        )

    def untag_resource(self, resource_arn: str, tag_keys: List[str]) -> None:
        self.tagger.untag_resource_using_names(resource_arn, tag_keys)

    def list_tags_for_resource(self, resource_arn: str) -> Dict[str, str]:
        return self.tagger.get_tag_dict_for_resource(resource_arn)

    def import_resources_to_draft_app_version(
        self,
        app_arn: str,
        eks_sources: List[Dict[str, Any]],
        source_arns: List[str],
        terraform_sources: List[Dict[str, str]],
    ) -> AppVersion:
        app = self.describe_app(app_arn)
        app_version = app.get_version("draft")

        app_version.eks_sources.extend(eks_sources)
        app_version.source_arns.extend(source_arns)
        app_version.terraform_sources.extend(terraform_sources)

        # Default AppComponent when importing data
        # AWS seems to create other components as well, based on the provided sources
        app_version.app_components.append(
            AppComponent(
                _id="appcommon",
                name="appcommon",
                _type="AWS::ResilienceHub::AppCommonAppComponent",
            )
        )
        return app_version

    def create_app_version_app_component(
        self, app_arn: str, name: str, _type: str
    ) -> AppComponent:
        app = self.describe_app(app_arn)
        app_version = app.get_version("draft")
        component = AppComponent(_id=name, name=name, _type=_type)
        app_version.app_components.append(component)
        return component

    def list_app_version_app_components(
        self, app_arn: str, app_version: str
    ) -> List[AppComponent]:
        app = self.describe_app(app_arn)
        return app.get_version(app_version).app_components

    def create_app_version_resource(
        self,
        app_arn: str,
        app_components: List[str],
        logical_resource_id: Dict[str, str],
        physical_resource_id: str,
        resource_type: str,
    ) -> Resource:
        app = self.describe_app(app_arn)
        app_version = app.get_version("draft")

        components = [c for c in app_version.app_components if c.id in app_components]

        resource = Resource(
            logical_resource_id=logical_resource_id,
            physical_resource_id=physical_resource_id,
            resource_type=resource_type,
            components=components,
        )
        app_version.resources.append(resource)
        return resource

    def list_app_version_resources(
        self, app_arn: str, app_version: str
    ) -> List[Resource]:
        app = self.describe_app(app_arn)
        return app.get_version(app_version).resources

    def list_app_versions(self, app_arn: str) -> List[AppVersion]:
        app = self.describe_app(app_arn)
        return app.app_versions

    def publish_app_version(self, app_arn: str, version_name: str) -> AppVersion:
        app = self.describe_app(app_arn)
        version = AppVersion(
            app_arn=app_arn, version_name=version_name, identifier=len(app.app_versions)
        )
        for old_version in app.app_versions:
            if old_version.app_version == "release":
                old_version.app_version = str(old_version.identifier)
        app.app_versions.append(version)
        return version


resiliencehub_backends = BackendDict(ResilienceHubBackend, "resiliencehub")
