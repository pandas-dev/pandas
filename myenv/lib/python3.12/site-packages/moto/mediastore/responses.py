import json

from moto.core.responses import BaseResponse

from .models import MediaStoreBackend, mediastore_backends


class MediaStoreResponse(BaseResponse):
    def __init__(self) -> None:
        super().__init__(service_name="mediastore")

    @property
    def mediastore_backend(self) -> MediaStoreBackend:
        return mediastore_backends[self.current_account][self.region]

    def create_container(self) -> str:
        name = self._get_param("ContainerName")
        tags = self._get_param("Tags")
        container = self.mediastore_backend.create_container(name=name, tags=tags)
        return json.dumps(dict(Container=container.to_dict()))

    def delete_container(self) -> str:
        name = self._get_param("ContainerName")
        self.mediastore_backend.delete_container(name=name)
        return "{}"

    def describe_container(self) -> str:
        name = self._get_param("ContainerName")
        container = self.mediastore_backend.describe_container(name=name)
        return json.dumps(dict(Container=container.to_dict()))

    def list_containers(self) -> str:
        containers = self.mediastore_backend.list_containers()
        return json.dumps(dict(dict(Containers=containers), NextToken=None))

    def list_tags_for_resource(self) -> str:
        name = self._get_param("Resource")
        tags = self.mediastore_backend.list_tags_for_resource(name)
        return json.dumps(dict(Tags=tags))

    def put_lifecycle_policy(self) -> str:
        container_name = self._get_param("ContainerName")
        lifecycle_policy = self._get_param("LifecyclePolicy")
        self.mediastore_backend.put_lifecycle_policy(
            container_name=container_name, lifecycle_policy=lifecycle_policy
        )
        return "{}"

    def get_lifecycle_policy(self) -> str:
        container_name = self._get_param("ContainerName")
        lifecycle_policy = self.mediastore_backend.get_lifecycle_policy(
            container_name=container_name
        )
        return json.dumps(dict(LifecyclePolicy=lifecycle_policy))

    def put_container_policy(self) -> str:
        container_name = self._get_param("ContainerName")
        policy = self._get_param("Policy")
        self.mediastore_backend.put_container_policy(
            container_name=container_name, policy=policy
        )
        return "{}"

    def get_container_policy(self) -> str:
        container_name = self._get_param("ContainerName")
        policy = self.mediastore_backend.get_container_policy(
            container_name=container_name
        )
        return json.dumps(dict(Policy=policy))

    def put_metric_policy(self) -> str:
        container_name = self._get_param("ContainerName")
        metric_policy = self._get_param("MetricPolicy")
        self.mediastore_backend.put_metric_policy(
            container_name=container_name, metric_policy=metric_policy
        )
        return json.dumps(metric_policy)

    def get_metric_policy(self) -> str:
        container_name = self._get_param("ContainerName")
        metric_policy = self.mediastore_backend.get_metric_policy(
            container_name=container_name
        )
        return json.dumps(dict(MetricPolicy=metric_policy))
