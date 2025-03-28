"""OpenSearchServiceServerlessBackend class with methods for supported APIs."""

import json
from typing import Any, Dict, List, Tuple

from moto.core.base_backend import BackendDict, BaseBackend
from moto.core.common_models import BaseModel
from moto.core.utils import unix_time
from moto.moto_api._internal import mock_random
from moto.utilities.tagging_service import TaggingService

from .exceptions import (
    ConflictException,
    ResourceNotFoundException,
    ValidationException,
)


class SecurityPolicy(BaseModel):
    def __init__(
        self,
        client_token: str,
        description: str,
        name: str,
        policy: str,
        type: str,
    ):
        self.client_token = client_token
        self.description = description
        self.name = name
        self.type = type
        self.created_date = int(unix_time() * 1000)
        # update policy # current date default
        self.last_modified_date = int(unix_time() * 1000)
        self.policy = json.loads(policy)
        self.policy_version = mock_random.get_random_string(20)
        if type == "encryption":
            self.resources = [
                res for rule in self.policy["Rules"] for res in rule["Resource"]
            ]
        else:
            self.resources = [
                res
                for p in self.policy
                for rule in p["Rules"]
                for res in rule["Resource"]
            ]

    def to_dict(self) -> Dict[str, Any]:
        dct = {
            "createdDate": self.created_date,
            "description": self.description,
            "lastModifiedDate": self.last_modified_date,
            "name": self.name,
            "policy": self.policy,
            "policyVersion": self.policy_version,
            "type": self.type,
        }
        return {k: v for k, v in dct.items() if v}

    def to_dict_list(self) -> Dict[str, Any]:
        dct = self.to_dict()
        dct.pop("policy")
        return {k: v for k, v in dct.items() if v}


class Collection(BaseModel):
    def __init__(
        self,
        client_token: str,
        description: str,
        name: str,
        standby_replicas: str,
        tags: List[Dict[str, str]],
        type: str,
        policy: Any,
        region: str,
        account_id: str,
    ):
        self.client_token = client_token
        self.description = description
        self.name = name
        self.standby_replicas = standby_replicas
        self.tags = tags
        self.type = type
        self.id = mock_random.get_random_string(length=20, lower_case=True)
        self.arn = f"arn:aws:aoss:{region}:{account_id}:collection/{self.id}"
        self.created_date = int(unix_time() * 1000)
        self.kms_key_arn = policy["KmsARN"]
        self.last_modified_date = int(unix_time() * 1000)
        self.status = "ACTIVE"
        self.collection_endpoint = f"https://{self.id}.{region}.aoss.amazonaws.com"
        self.dashboard_endpoint = (
            f"https://{self.id}.{region}.aoss.amazonaws.com/_dashboards"
        )

    def to_dict(self) -> Dict[str, Any]:
        dct = {
            "arn": self.arn,
            "createdDate": self.created_date,
            "description": self.description,
            "id": self.id,
            "kmsKeyArn": self.kms_key_arn,
            "lastModifiedDate": self.last_modified_date,
            "name": self.name,
            "standbyReplicas": self.standby_replicas,
            "status": self.status,
            "type": self.type,
        }
        return {k: v for k, v in dct.items() if v}

    def to_dict_list(self) -> Dict[str, Any]:
        dct = {"arn": self.arn, "id": self.id, "name": self.name, "status": self.status}
        return {k: v for k, v in dct.items() if v}

    def to_dict_batch(self) -> Dict[str, Any]:
        dct = self.to_dict()
        dct_options = {
            "collectionEndpoint": self.collection_endpoint,
            "dashboardEndpoint": self.dashboard_endpoint,
        }
        for key, value in dct_options.items():
            if value is not None:
                dct[key] = value
        return dct


class OSEndpoint(BaseModel):
    def __init__(
        self,
        client_token: str,
        name: str,
        security_group_ids: List[str],
        subnet_ids: List[str],
        vpc_id: str,
    ):
        self.client_token = client_token
        self.name = name
        self.security_group_ids = security_group_ids
        self.subnet_ids = subnet_ids
        self.vpc_id = vpc_id
        self.id = f"vpce-0{mock_random.get_random_string(length=16,lower_case=True)}"
        self.status = "ACTIVE"

    def to_dict(self) -> Dict[str, Any]:
        dct = {"id": self.id, "name": self.name, "status": self.status}
        return {k: v for k, v in dct.items() if v}


class OpenSearchServiceServerlessBackend(BaseBackend):
    """Implementation of OpenSearchServiceServerless APIs."""

    def __init__(self, region_name: str, account_id: str):
        super().__init__(region_name, account_id)

        self.collections: Dict[str, Collection] = dict()
        self.security_policies: Dict[str, SecurityPolicy] = dict()
        self.os_endpoints: Dict[str, OSEndpoint] = dict()
        self.tagger = TaggingService(
            tag_name="tags", key_name="key", value_name="value"
        )

    def create_security_policy(
        self, client_token: str, description: str, name: str, policy: str, type: str
    ) -> SecurityPolicy:
        if not client_token:
            client_token = mock_random.get_random_string(10)

        if (name, type) in list(
            (sp.name, sp.type) for sp in list(self.security_policies.values())
        ):
            raise ConflictException(
                msg=f"Policy with name {name} and type {type} already exists"
            )
        if type not in ["encryption", "network"]:
            raise ValidationException(
                msg=f"1 validation error detected: Value '{type}' at 'type' failed to satisfy constraint: Member must satisfy enum value set: [encryption, network]"
            )

        security_policy = SecurityPolicy(
            client_token=client_token,
            description=description,
            name=name,
            policy=policy,
            type=type,
        )
        self.security_policies[security_policy.client_token] = security_policy
        return security_policy

    def get_security_policy(self, name: str, type: str) -> SecurityPolicy:
        for sp in list(self.security_policies.values()):
            if sp.name == name and sp.type == type:
                return sp
        raise ResourceNotFoundException(
            msg=f"Policy with name {name} and type {type} is not found"
        )

    def list_security_policies(
        self, resource: List[str], type: str
    ) -> List[SecurityPolicy]:
        """
        Pagination is not yet implemented
        """
        security_policy_summaries = []
        if resource:
            for res in resource:
                security_policy_summaries.extend(
                    [
                        sp
                        for sp in list(self.security_policies.values())
                        if res in sp.resources and type == sp.type
                    ]
                )
        else:
            security_policy_summaries = [
                sp for sp in list(self.security_policies.values()) if sp.type == type
            ]
        return security_policy_summaries

    def update_security_policy(
        self,
        client_token: str,
        description: str,
        name: str,
        policy: str,
        policy_version: str,
        type: str,
    ) -> SecurityPolicy:
        if not client_token:
            client_token = mock_random.get_random_string(10)

        for sp in list(self.security_policies.values()):
            if sp.name == name and sp.type == type:
                if sp.policy_version == policy_version:
                    last_modified_date = sp.last_modified_date
                    if sp.policy != json.loads(policy):
                        last_modified_date = int(unix_time() * 1000)
                        # Updating policy version
                        policy_version = mock_random.get_random_string(20)

                    sp.client_token = client_token
                    sp.description = description
                    sp.name = name
                    sp.policy = json.loads(policy)
                    sp.last_modified_date = last_modified_date
                    sp.policy_version = policy_version
                    return sp
                else:
                    raise ValidationException(
                        msg="Policy version specified in the request refers to an older version and policy has since changed"
                    )

        raise ResourceNotFoundException(
            msg=f"Policy with name {name} and type {type} is not found"
        )

    def create_collection(
        self,
        client_token: str,
        description: str,
        name: str,
        standby_replicas: str,
        tags: List[Dict[str, str]],
        type: str,
    ) -> Collection:
        policy = ""
        if not client_token:
            client_token = mock_random.get_random_string(10)

        for sp in list(self.security_policies.values()):
            if f"collection/{name}" in sp.resources:
                policy = sp.policy
        if not policy:
            raise ValidationException(
                msg=f"No matching security policy of encryption type found for collection name: {name}. Please create security policy of encryption type for this collection."
            )

        collection = Collection(
            client_token=client_token,
            description=description,
            name=name,
            standby_replicas=standby_replicas,
            tags=tags,
            type=type,
            policy=policy,
            region=self.region_name,
            account_id=self.account_id,
        )
        self.collections[collection.id] = collection
        self.tag_resource(collection.arn, tags)
        return collection

    def list_collections(self, collection_filters: Dict[str, str]) -> List[Collection]:
        """
        Pagination is not yet implemented
        """
        collection_summaries = []
        if (collection_filters) and ("name" in collection_filters):
            collection_summaries = [
                collection
                for collection in list(self.collections.values())
                if collection.name == collection_filters["name"]
            ]
        else:
            collection_summaries = [
                collection for collection in list(self.collections.values())
            ]
        return collection_summaries

    def create_vpc_endpoint(
        self,
        client_token: str,
        name: str,
        security_group_ids: List[str],
        subnet_ids: List[str],
        vpc_id: str,
    ) -> OSEndpoint:
        if not client_token:
            client_token = mock_random.get_random_string(10)

        # Only 1 endpoint should exists under each VPC
        if vpc_id in list(ose.vpc_id for ose in list(self.os_endpoints.values())):
            raise ConflictException(
                msg=f"Failed to create a VpcEndpoint {name} for AccountId {self.account_id} :: There is already a VpcEndpoint exist under VpcId {vpc_id}"
            )

        os_endpoint = OSEndpoint(
            client_token=client_token,
            name=name,
            security_group_ids=security_group_ids,
            subnet_ids=subnet_ids,
            vpc_id=vpc_id,
        )
        self.os_endpoints[os_endpoint.client_token] = os_endpoint

        return os_endpoint

    def delete_collection(self, client_token: str, id: str) -> Collection:
        if not client_token:
            client_token = mock_random.get_random_string(10)

        if id in self.collections:
            self.collections[id].status = "DELETING"
            return self.collections.pop(id)
        raise ResourceNotFoundException(f"Collection with ID {id} cannot be found.")

    def tag_resource(self, resource_arn: str, tags: List[Dict[str, str]]) -> None:
        self.tagger.tag_resource(resource_arn, tags)

    def untag_resource(self, resource_arn: str, tag_keys: List[str]) -> None:
        self.tagger.untag_resource_using_names(resource_arn, tag_keys)

    def list_tags_for_resource(self, resource_arn: str) -> List[Dict[str, str]]:
        return self.tagger.list_tags_for_resource(resource_arn)["tags"]

    def batch_get_collection(
        self, ids: List[str], names: List[str]
    ) -> Tuple[List[Any], List[Dict[str, str]]]:
        collection_details = []
        collection_error_details = []
        collection_error_detail = {
            "errorCode": "NOT_FOUND",
            "errorMessage": "The specified Collection is not found.",
        }
        if ids and names:
            raise ValidationException(
                msg="You need to provide IDs or names. You can't provide both IDs and names in the same request"
            )
        if ids:
            for i in ids:
                if i in self.collections:
                    collection_details.append(self.collections[i].to_dict_batch())
                else:
                    collection_error_detail["id"] = i
                    collection_error_details.append(collection_error_detail)

        if names:
            for n in names:
                for collection in self.collections.values():
                    if collection.name == n:
                        collection_details.append(collection.to_dict_batch())
                    else:
                        collection_error_detail["name"] = n
                        collection_error_details.append(collection_error_detail)
        return collection_details, collection_error_details


opensearchserverless_backends = BackendDict(
    OpenSearchServiceServerlessBackend, "opensearchserverless"
)
