"""S3TablesBackend class with methods for supported APIs."""

import datetime
import re
from hashlib import md5
from typing import Dict, List, Literal, Optional, Union

from moto.core.base_backend import BackendDict, BaseBackend
from moto.moto_api._internal import mock_random as random
from moto.s3.models import FakeBucket
from moto.s3tables.exceptions import (
    ConflictException,
    DestinationNamespaceDoesNotExist,
    InvalidContinuationToken,
    InvalidMetadataLocation,
    InvalidNamespaceName,
    InvalidTableBucketName,
    InvalidTableName,
    NotFoundException,
    NothingToRename,
    TableAlreadyExists,
    TableDoesNotExist,
    VersionTokenMismatch,
)
from moto.utilities.paginator import paginate
from moto.utilities.utils import get_partition

# https://docs.aws.amazon.com/AmazonS3/latest/userguide/s3-tables-buckets-naming.html
TABLE_BUCKET_NAME_PATTERN = re.compile(r"[a-z0-9_-]{3,63}")
TABLE_BUCKET_NAME_RESERVED_PREFIXES = ("xn--", "sthree-", "amzn-s3-demo")
TABLE_BUCKET_NAME_RESERVED_SUFFIXES = ("-s3alias", "--ol-s3", "--x-s3")
NAMESPACE_NAME_PATTERN = re.compile(r"[0-9a-z_]*")
TABLE_NAME_PATTERN = re.compile(r"[0-9a-z_]*")


def _validate_table_bucket_name(name: str) -> None:
    if (
        not TABLE_BUCKET_NAME_PATTERN.match(name)
        or any(
            name.startswith(prefix) for prefix in TABLE_BUCKET_NAME_RESERVED_PREFIXES
        )
        or any(name.endswith(suffix) for suffix in TABLE_BUCKET_NAME_RESERVED_SUFFIXES)
    ):
        raise InvalidTableBucketName()


def _validate_namespace_name(name: str) -> None:
    if not NAMESPACE_NAME_PATTERN.match(name):
        raise InvalidNamespaceName()


def _validate_table_name(name: str) -> None:
    if not TABLE_NAME_PATTERN.match(name):
        raise InvalidTableName(name)


PAGINATION_MODEL = {
    "list_table_buckets": {
        "input_token": "continuation_token",
        "limit_key": "max_buckets",
        "limit_default": 1000,
        "unique_attribute": ["arn"],
        "fail_on_invalid_token": InvalidContinuationToken,
    },
    "list_namespaces": {
        "input_token": "continuation_token",
        "limit_key": "max_namespaces",
        "limit_default": 1000,
        "unique_attribute": ["name"],
        "fail_on_invalid_token": InvalidContinuationToken,
    },
    "list_tables": {
        "input_token": "continuation_token",
        "limit_key": "max_tables",
        "limit_default": 1000,
        "unique_attribute": ["arn"],
        "fail_on_invalid_token": InvalidContinuationToken,
    },
}


class Table:
    def __init__(
        self,
        name: str,
        account_id: str,
        created_by: str,
        format: Literal["ICEBERG"],
        namespace: str,
        table_bucket_arn: str,
        type: Union[Literal["customer"], Literal["aws"]],
        managed_by_service: bool,
        partition: str,
    ):
        _validate_table_name(name)
        self.name = name
        self.account_id = account_id
        self.partition = partition
        self.created_by = created_by
        self.format = format
        self.type = type
        self.version_token = self._generate_version_token()
        self.creation_date = datetime.datetime.now(tz=datetime.timezone.utc)
        self.last_modified = self.creation_date
        self.modified_by: Optional[str] = None
        self.namespace = namespace
        self.table_bucket_arn = table_bucket_arn
        self.region_name = table_bucket_arn.split(":")[3]
        self.managed_by_service = managed_by_service
        self.metadata_location: Optional[str] = None
        self._bucket = self._create_underlying_bucket()
        self.warehouse_location: str = f"s3://{self._bucket.name}"

    @property
    def arn(self) -> str:
        return f"{self.table_bucket_arn}/table/{self.name}"

    def was_modified(self, by: str) -> None:
        self.last_modified = datetime.datetime.now(tz=datetime.timezone.utc)
        self.modified_by = by

    def _generate_version_token(self) -> str:
        return md5(random.uuid4().bytes).hexdigest()[:20]

    def _create_underlying_bucket(self) -> FakeBucket:
        from moto.s3.models import s3_backends

        bucket = s3_backends[self.account_id][
            self.partition
        ].create_table_storage_bucket(region_name=self.region_name)
        return bucket

    def update_metadata_location(
        self, metadata_location: str, version_token: str
    ) -> None:
        if not metadata_location.startswith(self.warehouse_location):
            raise InvalidMetadataLocation()
        if not self.version_token == version_token:
            raise VersionTokenMismatch()

        self.metadata_location = metadata_location
        self.version_token = self._generate_version_token()

    def rename(self, new_name: str, by: str) -> None:
        _validate_table_name(new_name)
        self.name = new_name
        self.was_modified(by)


class Namespace:
    def __init__(self, name: str, account_id: str, created_by: str):
        _validate_namespace_name(name)
        self.name = name
        self.account_id = account_id
        self.created_by = created_by
        self.creation_date = datetime.datetime.now(tz=datetime.timezone.utc)
        self.tables: Dict[str, Table] = {}


class FakeTableBucket:
    def __init__(self, name: str, account_id: str, region_name: str):
        self.name = name
        self.account_id = account_id
        self.region_name = region_name
        self.partition = get_partition(region_name)
        self.creation_date = datetime.datetime.now(tz=datetime.timezone.utc)
        self.namespaces: Dict[str, Namespace] = {}

    @property
    def arn(self) -> str:
        return f"arn:{self.partition}:s3tables:{self.region_name}:{self.account_id}:bucket/{self.name}"


class S3TablesBackend(BaseBackend):
    """Implementation of S3Tables APIs."""

    def __init__(self, region_name: str, account_id: str) -> None:
        super().__init__(region_name, account_id)
        self.table_buckets: Dict[str, FakeTableBucket] = {}

    def create_table_bucket(self, name: str) -> FakeTableBucket:
        _validate_table_bucket_name(name)
        new_table_bucket = FakeTableBucket(
            name=name, account_id=self.account_id, region_name=self.region_name
        )
        if new_table_bucket.arn in self.table_buckets:
            raise ConflictException(
                "The bucket that you tried to create already exists, and you own it."
            )
        self.table_buckets[new_table_bucket.arn] = new_table_bucket

        return new_table_bucket

    @paginate(pagination_model=PAGINATION_MODEL)
    def list_table_buckets(
        self,
        prefix: Optional[str] = None,
    ) -> List[FakeTableBucket]:
        all_buckets = list(
            bucket
            for bucket in self.table_buckets.values()
            if (prefix is None or bucket.name.startswith(prefix))
        )

        return all_buckets

    def get_table_bucket(self, table_bucket_arn: str) -> FakeTableBucket:
        bucket = self.table_buckets.get(table_bucket_arn)
        if not bucket:
            raise NotFoundException("The specified bucket does not exist.")
        return bucket

    def delete_table_bucket(self, table_bucket_arn: str) -> None:
        # make sure table exists first
        self.get_table_bucket(table_bucket_arn)
        self.table_buckets.pop(table_bucket_arn)

    def create_namespace(self, table_bucket_arn: str, namespace: str) -> Namespace:
        bucket = self.table_buckets.get(table_bucket_arn)

        if not bucket:
            raise NotFoundException(
                "The request was rejected because the specified resource could not be found."
            )

        if namespace in bucket.namespaces:
            raise ConflictException(
                "A namespace with an identical name already exists in the bucket."
            )

        ns = Namespace(
            namespace, account_id=self.account_id, created_by=self.account_id
        )
        bucket.namespaces[ns.name] = ns
        return ns

    @paginate(pagination_model=PAGINATION_MODEL)
    def list_namespaces(
        self,
        table_bucket_arn: str,
        prefix: Optional[str] = None,
    ) -> List[Namespace]:
        bucket = self.get_table_bucket(table_bucket_arn)

        all_namespaces = list(
            ns
            for ns in bucket.namespaces.values()
            if (prefix is None or ns.name.startswith(prefix))
        )

        return all_namespaces

    def get_namespace(self, table_bucket_arn: str, namespace: str) -> Namespace:
        bucket = self.table_buckets.get(table_bucket_arn)
        if bucket and namespace in bucket.namespaces:
            return bucket.namespaces[namespace]

        raise NotFoundException(
            "The request was rejected because the specified resource could not be found."
        )

    def delete_namespace(self, table_bucket_arn: str, namespace: str) -> None:
        bucket = self.table_buckets.get(table_bucket_arn)
        if bucket and namespace in bucket.namespaces:
            bucket.namespaces.pop(namespace)
            return

        raise NotFoundException(
            "The request was rejected because the specified resource could not be found."
        )

    def create_table(
        self,
        table_bucket_arn: str,
        namespace: str,
        name: str,
        format: Literal["ICEBERG"],
    ) -> Table:
        bucket = self.table_buckets.get(table_bucket_arn)
        if not bucket:
            raise NotFoundException("The specified bucket does not exist.")

        if namespace not in bucket.namespaces:
            raise NotFoundException("The specified namespace does not exist.")

        ns = bucket.namespaces[namespace]
        if name in ns.tables:
            TableAlreadyExists()
        table = Table(
            name=name,
            account_id=self.account_id,
            created_by=self.account_id,
            format=format,
            namespace=namespace,
            table_bucket_arn=table_bucket_arn,
            type="customer",
            managed_by_service=False,
            partition=self.partition,
        )
        ns.tables[table.name] = table
        return table

    def get_table(self, table_bucket_arn: str, namespace: str, name: str) -> Table:
        bucket = self.table_buckets.get(table_bucket_arn)
        if bucket and namespace in bucket.namespaces:
            if name in bucket.namespaces[namespace].tables:
                return bucket.namespaces[namespace].tables[name]
        raise TableDoesNotExist()

    @paginate(pagination_model=PAGINATION_MODEL)
    def list_tables(
        self,
        table_bucket_arn: str,
        namespace: Optional[str] = None,
        prefix: Optional[str] = None,
    ) -> List[Table]:
        bucket = self.table_buckets.get(table_bucket_arn)
        if not bucket or (namespace and namespace not in bucket.namespaces):
            raise NotFoundException(
                "The request was rejected because the specified resource could not be found."
            )

        if namespace:
            all_tables = list(
                table
                for table in bucket.namespaces[namespace].tables.values()
                if (prefix is None or table.name.startswith(prefix))
            )
        else:
            all_tables = list(
                table
                for ns in bucket.namespaces.values()
                for table in ns.tables.values()
                if (prefix is None or table.name.startswith(prefix))
            )

        return all_tables

    def delete_table(
        self,
        table_bucket_arn: str,
        namespace: str,
        name: str,
        version_token: Optional[str] = None,
    ) -> None:
        bucket = self.table_buckets.get(table_bucket_arn)
        if (
            not bucket
            or namespace not in bucket.namespaces
            or (name not in bucket.namespaces[namespace].tables)
        ):
            raise TableDoesNotExist()

        ns = bucket.namespaces[namespace]
        table = ns.tables[name]
        if version_token and not version_token == table.version_token:
            raise VersionTokenMismatch()
        from moto.s3.models import s3_backends

        s3_backends[self.account_id][self.partition].delete_table_storage_bucket(
            table._bucket.name
        )

        ns.tables.pop(name)

    def update_table_metadata_location(
        self,
        table_bucket_arn: str,
        namespace: str,
        name: str,
        version_token: str,
        metadata_location: str,
    ) -> Table:
        bucket = self.table_buckets.get(table_bucket_arn)
        if bucket and namespace in bucket.namespaces:
            if name in bucket.namespaces[namespace].tables:
                table = bucket.namespaces[namespace].tables[name]
                table.update_metadata_location(
                    metadata_location=metadata_location, version_token=version_token
                )
                return table

        raise TableDoesNotExist()

    def rename_table(
        self,
        table_bucket_arn: str,
        namespace: str,
        name: str,
        new_namespace_name: Optional[str] = None,
        new_name: Optional[str] = None,
        version_token: Optional[str] = None,
    ) -> None:
        if not new_namespace_name and not new_name:
            raise NothingToRename()
        destination_namespace = new_namespace_name if new_namespace_name else namespace
        destination_name = new_name if new_name else name
        _validate_table_name(destination_name)

        bucket = self.table_buckets.get(table_bucket_arn)
        if not bucket or destination_namespace not in bucket.namespaces:
            raise DestinationNamespaceDoesNotExist()
        if namespace not in bucket.namespaces or (
            name not in bucket.namespaces[namespace].tables
        ):
            raise TableDoesNotExist()
        table = bucket.namespaces[namespace].tables[name]
        if version_token and not version_token == table.version_token:
            raise VersionTokenMismatch()

        if destination_name in bucket.namespaces[destination_namespace].tables:
            raise TableAlreadyExists()

        table = bucket.namespaces[namespace].tables.pop(name)
        table.rename(new_name=destination_name, by=self.account_id)
        bucket.namespaces[destination_namespace].tables[destination_name] = table


s3tables_backends = BackendDict(
    S3TablesBackend,
    "s3tables",
    additional_regions=["us-east-1", "us-east-2", "us-west-2"],
)
