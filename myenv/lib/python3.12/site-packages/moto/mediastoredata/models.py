import hashlib
from collections import OrderedDict
from typing import Any, Dict, List

from moto.core.base_backend import BackendDict, BaseBackend
from moto.core.common_models import BaseModel

from .exceptions import ClientError


class Object(BaseModel):
    def __init__(
        self, path: str, body: str, etag: str, storage_class: str = "TEMPORAL"
    ):
        self.path = path
        self.body = body
        self.content_sha256 = hashlib.sha256(body.encode("utf-8")).hexdigest()
        self.etag = etag
        self.storage_class = storage_class

    def to_dict(self) -> Dict[str, Any]:
        return {
            "ETag": self.etag,
            "Name": self.path,
            "Type": "FILE",
            "ContentLength": 123,
            "StorageClass": self.storage_class,
            "Path": self.path,
            "ContentSHA256": self.content_sha256,
        }


class MediaStoreDataBackend(BaseBackend):
    def __init__(self, region_name: str, account_id: str):
        super().__init__(region_name, account_id)
        self._objects: Dict[str, Object] = OrderedDict()

    def put_object(
        self, body: str, path: str, storage_class: str = "TEMPORAL"
    ) -> Object:
        """
        The following parameters are not yet implemented: ContentType, CacheControl, UploadAvailability
        """
        new_object = Object(
            path=path, body=body, etag="etag", storage_class=storage_class
        )
        self._objects[path] = new_object
        return new_object

    def delete_object(self, path: str) -> None:
        if path not in self._objects:
            raise ClientError(
                "ObjectNotFoundException", f"Object with id={path} not found"
            )
        del self._objects[path]

    def get_object(self, path: str) -> Object:
        """
        The Range-parameter is not yet supported.
        """
        objects_found = [item for item in self._objects.values() if item.path == path]
        if len(objects_found) == 0:
            raise ClientError(
                "ObjectNotFoundException", f"Object with id={path} not found"
            )
        return objects_found[0]

    def list_items(self) -> List[Dict[str, Any]]:
        """
        The Path- and MaxResults-parameters are not yet supported.
        """
        return [c.to_dict() for c in self._objects.values()]


mediastoredata_backends = BackendDict(MediaStoreDataBackend, "mediastore-data")
