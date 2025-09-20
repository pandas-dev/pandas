import functools
from typing import TYPE_CHECKING, Any, Tuple

import requests.adapters

from moto import settings

if TYPE_CHECKING:
    from docker import DockerClient


_orig_adapter_send = requests.adapters.HTTPAdapter.send


class DockerModel:
    def __init__(self) -> None:
        self.__docker_client = None

    @property
    def docker_client(self) -> "DockerClient":  # type: ignore
        if self.__docker_client is None:
            # We should only initiate the Docker Client at runtime.
            # The docker.from_env() call will fall if Docker is not running
            import docker

            self.__docker_client = docker.from_env()

            # Unfortunately mocking replaces this method w/o fallback enabled, so we
            # need to replace it if we detect it's been mocked
            if requests.adapters.HTTPAdapter.send != _orig_adapter_send:
                _orig_get_adapter = self.docker_client.api.get_adapter

                def replace_adapter_send(*args: Any, **kwargs: Any) -> Any:
                    adapter = _orig_get_adapter(*args, **kwargs)

                    if isinstance(adapter, requests.adapters.HTTPAdapter):
                        adapter.send = functools.partial(_orig_adapter_send, adapter)  # type: ignore
                    return adapter

                self.docker_client.api.get_adapter = replace_adapter_send
        return self.__docker_client

    def ensure_image_exists(self, name: str) -> None:
        full_name = ":".join(parse_image_ref(name))
        try:
            self.docker_client.images.get(full_name)
        except:  # noqa: E722 Do not use bare except
            self.docker_client.images.pull(full_name)


def parse_image_ref(image_name: str) -> Tuple[str, str]:
    # podman does not support short container image name out of box - try to make a full name
    # See ParseDockerRef() in https://github.com/distribution/distribution/blob/main/reference/normalize.go
    parts = image_name.split("/")
    if len(parts) == 1 or (
        "." not in parts[0] and ":" not in parts[0] and parts[0] != "localhost"
    ):
        domain = settings.DEFAULT_CONTAINER_REGISTRY
        remainder = parts
    else:
        domain = parts[0]
        remainder = parts[1:]
    # Special handling for docker.io
    # https://github.com/containers/image/blob/master/docs/containers-registries.conf.5.md#normalization-of-dockerio-references
    if domain == "docker.io" and len(remainder) == 1:
        remainder = ["library"] + remainder
    if ":" in remainder[-1]:
        remainder[-1], image_tag = remainder[-1].split(":", 1)
    else:
        image_tag = "latest"
    image_repository = "/".join([domain] + remainder)
    return image_repository, image_tag
