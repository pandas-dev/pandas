from collections import OrderedDict
from datetime import date
from typing import Any, Dict, List, Optional

from moto.core.base_backend import BackendDict, BaseBackend
from moto.core.common_models import BaseModel
from moto.utilities.utils import get_partition

from .exceptions import (
    ContainerNotFoundException,
    PolicyNotFoundException,
    ResourceNotFoundException,
)


class Container(BaseModel):
    def __init__(self, **kwargs: Any):
        self.arn = kwargs.get("arn")
        self.name = kwargs.get("name")
        self.endpoint = kwargs.get("endpoint")
        self.status = kwargs.get("status")
        self.creation_time = kwargs.get("creation_time")
        self.lifecycle_policy: Optional[str] = None
        self.policy: Optional[str] = None
        self.metric_policy: Optional[str] = None
        self.tags = kwargs.get("tags")

    def to_dict(self, exclude: Optional[List[str]] = None) -> Dict[str, Any]:
        data = {
            "ARN": self.arn,
            "Name": self.name,
            "Endpoint": self.endpoint,
            "Status": self.status,
            "CreationTime": self.creation_time,
            "Tags": self.tags,
        }
        if exclude:
            for key in exclude:
                del data[key]
        return data


class MediaStoreBackend(BaseBackend):
    def __init__(self, region_name: str, account_id: str):
        super().__init__(region_name, account_id)
        self._containers: Dict[str, Container] = OrderedDict()

    def create_container(self, name: str, tags: Dict[str, str]) -> Container:
        arn = f"arn:{get_partition(self.region_name)}:mediastore:container:{name}"
        container = Container(
            arn=arn,
            name=name,
            endpoint=f"/{name}",
            status="CREATING",
            creation_time=date.today().strftime("%m/%d/%Y, %H:%M:%S"),
            tags=tags,
        )
        self._containers[name] = container
        return container

    def delete_container(self, name: str) -> None:
        if name not in self._containers:
            raise ContainerNotFoundException()
        del self._containers[name]

    def describe_container(self, name: str) -> Container:
        if name not in self._containers:
            raise ResourceNotFoundException()
        container = self._containers[name]
        container.status = "ACTIVE"
        return container

    def list_containers(self) -> List[Dict[str, Any]]:
        """
        Pagination is not yet implemented
        """
        return [c.to_dict() for c in self._containers.values()]

    def list_tags_for_resource(self, name: str) -> Optional[Dict[str, str]]:
        if name not in self._containers:
            raise ContainerNotFoundException()
        tags = self._containers[name].tags
        return tags

    def put_lifecycle_policy(self, container_name: str, lifecycle_policy: str) -> None:
        if container_name not in self._containers:
            raise ResourceNotFoundException()
        self._containers[container_name].lifecycle_policy = lifecycle_policy

    def get_lifecycle_policy(self, container_name: str) -> str:
        if container_name not in self._containers:
            raise ResourceNotFoundException()
        lifecycle_policy = self._containers[container_name].lifecycle_policy
        if not lifecycle_policy:
            raise PolicyNotFoundException()
        return lifecycle_policy

    def put_container_policy(self, container_name: str, policy: str) -> None:
        if container_name not in self._containers:
            raise ResourceNotFoundException()
        self._containers[container_name].policy = policy

    def get_container_policy(self, container_name: str) -> str:
        if container_name not in self._containers:
            raise ResourceNotFoundException()
        policy = self._containers[container_name].policy
        if not policy:
            raise PolicyNotFoundException()
        return policy

    def put_metric_policy(self, container_name: str, metric_policy: str) -> None:
        if container_name not in self._containers:
            raise ResourceNotFoundException()
        self._containers[container_name].metric_policy = metric_policy

    def get_metric_policy(self, container_name: str) -> str:
        if container_name not in self._containers:
            raise ResourceNotFoundException()
        metric_policy = self._containers[container_name].metric_policy
        if not metric_policy:
            raise PolicyNotFoundException()
        return metric_policy


mediastore_backends = BackendDict(MediaStoreBackend, "mediastore")
