"""S3VectorsBackend class with methods for supported APIs."""

from collections.abc import Iterable
from typing import Any, Literal, Optional, TypedDict

from moto.core.base_backend import BackendDict, BaseBackend
from moto.core.common_models import BaseModel
from moto.utilities.arns import parse_arn
from moto.utilities.utils import PARTITION_NAMES

from .exceptions import (
    IndexNotFound,
    VectorBucketAlreadyExists,
    VectorBucketNotEmpty,
    VectorBucketNotFound,
    VectorBucketPolicyNotFound,
    VectorWrongDimension,
)
from .utils import create_vector_bucket_arn


class VectorData(TypedDict):
    float32: list[float]


class VectorType(TypedDict, total=False):
    key: str
    data: VectorData
    metadata: Any


class Vector(BaseModel):
    def __init__(self, key: str, data: VectorData, metadata: Any):
        self.key = key
        self.data = data
        self.metadata = metadata

    def to_dict(self, return_data: bool, return_metadata: bool) -> VectorType:
        return (
            {"key": self.key}  # type: ignore[return-value]
            | ({"data": self.data} if return_data else {})
            | ({"metadata": self.metadata} if return_metadata else {})
        )

    @staticmethod
    def from_dict(dct: VectorType) -> "Vector":
        return Vector(
            key=dct["key"], data=dct["data"], metadata=dct.get("metadata", {})
        )


class Index(BaseModel):
    def __init__(
        self,
        bucket: "VectorBucket",
        name: str,
        dimension: int,
        data_type: Literal["float32"],
        distance_metric: str,
    ):
        self.vectorBucketName = bucket.vector_bucket_name
        self.index_name = name
        self.index_arn = f"{bucket.vector_bucket_arn}/index/{name}"
        self.dimension = dimension
        self.data_type = data_type
        self.distance_metric = distance_metric

        self._bucket = bucket
        self.vectors: dict[str, Vector] = {}


class VectorBucket(BaseModel):
    def __init__(
        self,
        arn: str,
        name: str,
        encryption_configuration: dict[str, str],
    ):
        self.vector_bucket_name = name
        self.vector_bucket_arn = arn
        self.encryption_configuration = encryption_configuration or {
            "sseType": "AES256"
        }

        self.indexes: dict[str, Index] = {}
        self.policy: Optional[str] = None


class S3VectorsBackend(BaseBackend):
    """Implementation of S3Vectors APIs."""

    def __init__(self, region_name: str, account_id: str):
        super().__init__(region_name, account_id)
        self.vector_buckets: dict[str, VectorBucket] = {}

    def create_vector_bucket(
        self,
        region: str,
        vector_bucket_name: str,
        encryption_configuration: dict[str, str],
    ) -> None:
        vector_bucket_arn = create_vector_bucket_arn(
            self.account_id, region, name=vector_bucket_name
        )
        if vector_bucket_arn in self.vector_buckets:
            raise VectorBucketAlreadyExists
        vector_bucket = VectorBucket(
            arn=vector_bucket_arn,
            name=vector_bucket_name,
            encryption_configuration=encryption_configuration,
        )
        self.vector_buckets[vector_bucket.vector_bucket_arn] = vector_bucket

    def get_vector_bucket(
        self,
        vector_bucket_name: Optional[str] = None,
        vector_bucket_arn: Optional[str] = None,
    ) -> VectorBucket:
        if vector_bucket_name:
            for vector_bucket in self.vector_buckets.values():
                if vector_bucket.vector_bucket_name == vector_bucket_name:
                    return vector_bucket
        if vector_bucket_arn and (bucket := self.vector_buckets.get(vector_bucket_arn)):
            return bucket
        raise VectorBucketNotFound

    def delete_vector_bucket(self, vector_bucket_name: str) -> None:
        if vector_bucket_name:
            bucket = self.get_vector_bucket(vector_bucket_name=vector_bucket_name)
            if bucket.indexes:
                raise VectorBucketNotEmpty
            self.vector_buckets.pop(bucket.vector_bucket_arn, None)

    def list_vector_buckets(self, prefix: Optional[str]) -> list[VectorBucket]:
        return [
            bucket
            for bucket in self.vector_buckets.values()
            if not prefix or bucket.vector_bucket_name.startswith(prefix)
        ]

    def create_index(
        self,
        vector_bucket_name: str,
        vector_bucket_arn: str,
        index_name: str,
        data_type: Literal["float32"],
        dimension: int,
        distance_metric: str,
    ) -> None:
        bucket = self.get_vector_bucket(
            vector_bucket_name=vector_bucket_name, vector_bucket_arn=vector_bucket_arn
        )
        index = Index(
            bucket=bucket,
            name=index_name,
            data_type=data_type,
            dimension=dimension,
            distance_metric=distance_metric,
        )
        bucket.indexes[index.index_arn] = index

    def delete_index(
        self, vector_bucket_name: str, index_name: str, index_arn: str
    ) -> None:
        index = self.get_index(vector_bucket_name, index_name, index_arn)
        index._bucket.indexes.pop(index.index_arn)

    def get_index(
        self, vector_bucket_name: str, index_name: str, index_arn: str
    ) -> Index:
        if index_arn:
            vector_bucket_name, _, index_name = parse_arn(index_arn).resource_id.split(
                "/"
            )
        try:
            bucket = self.get_vector_bucket(
                vector_bucket_name=vector_bucket_name, vector_bucket_arn=None
            )
            for index in bucket.indexes.values():
                if index.index_name == index_name:
                    return index
        except VectorBucketNotFound:
            pass
        raise IndexNotFound

    def list_indexes(
        self, vector_bucket_name: str, vector_bucket_arn: str
    ) -> list[Index]:
        """Pagination is not yet implemented. The prefix-parameter is also not yet implemented."""
        bucket = self.get_vector_bucket(
            vector_bucket_name, vector_bucket_arn=vector_bucket_arn
        )
        return list(bucket.indexes.values())

    def delete_vector_bucket_policy(
        self, vector_bucket_name: str, vector_bucket_arn: str
    ) -> None:
        bucket = self.get_vector_bucket(
            vector_bucket_name, vector_bucket_arn=vector_bucket_arn
        )
        bucket.policy = None

    def get_vector_bucket_policy(
        self, vector_bucket_name: str, vector_bucket_arn: str
    ) -> str:
        bucket = self.get_vector_bucket(
            vector_bucket_name, vector_bucket_arn=vector_bucket_arn
        )
        if not bucket.policy:
            raise VectorBucketPolicyNotFound
        return bucket.policy

    def put_vector_bucket_policy(
        self, vector_bucket_name: str, vector_bucket_arn: str, policy: str
    ) -> None:
        bucket = self.get_vector_bucket(
            vector_bucket_name, vector_bucket_arn=vector_bucket_arn
        )
        bucket.policy = policy

    def put_vectors(
        self,
        vector_bucket_name: str,
        index_name: str,
        index_arn: str,
        vectors: list[VectorType],
    ) -> None:
        index = self.get_index(
            vector_bucket_name, index_name=index_name, index_arn=index_arn
        )

        for vector in vectors:
            provided = len(vector["data"][index.data_type])  # type: ignore[literal-required]
            if provided != index.dimension:
                raise VectorWrongDimension(
                    key=vector["key"], actual=index.dimension, provided=provided
                )

        for vector in vectors:
            index.vectors[vector["key"]] = Vector.from_dict(vector)

    def get_vectors(
        self, vector_bucket_name: str, index_name: str, index_arn: str, keys: list[str]
    ) -> list[Vector]:
        index = self.get_index(
            vector_bucket_name, index_name=index_name, index_arn=index_arn
        )
        return [index.vectors[name] for name in index.vectors if name in keys]

    def list_vectors(
        self, vector_bucket_name: str, index_name: str, index_arn: str
    ) -> Iterable[Vector]:
        """
        Pagination is not yet implemented
        Segmentation is not yet implemented
        """
        index = self.get_index(
            vector_bucket_name, index_name=index_name, index_arn=index_arn
        )
        return index.vectors.values()

    def delete_vectors(
        self, vector_bucket_name: str, index_name: str, index_arn: str, keys: list[str]
    ) -> None:
        index = self.get_index(
            vector_bucket_name, index_name=index_name, index_arn=index_arn
        )
        for key in keys:
            index.vectors.pop(key, None)


s3vectors_backends = BackendDict(
    S3VectorsBackend,
    "s3vectors",
    use_boto3_regions=False,
    additional_regions=PARTITION_NAMES,
)
