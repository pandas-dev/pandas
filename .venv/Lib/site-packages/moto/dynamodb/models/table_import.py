from threading import Thread
from typing import TYPE_CHECKING, Any, Dict, List, Optional
from uuid import uuid4

from moto.core.utils import gzip_decompress
from moto.utilities.utils import get_partition

if TYPE_CHECKING:
    from moto.dynamodb.models import DynamoDBBackend
    from moto.dynamodb.models.table import Table
    from moto.s3.models import FakeBucket, S3Backend


class TableImport(Thread):
    def __init__(
        self,
        account_id: str,
        s3_source: Dict[str, str],
        region_name: str,
        table_name: str,
        billing_mode: str,
        throughput: Optional[Dict[str, int]],
        key_schema: List[Dict[str, str]],
        global_indexes: Optional[List[Dict[str, Any]]],
        attrs: List[Dict[str, str]],
        compression_type: Optional[str],
    ):
        super().__init__()
        self.partition = get_partition(region_name)
        self.arn = f"arn:{self.partition}:dynamodb:{region_name}:{account_id}:table/{table_name}/import/{str(uuid4()).replace('-', '')}"
        self.status = "IN_PROGRESS"
        self.account_id = account_id
        self.s3_source = s3_source
        self.region_name = region_name

        self.table_name = table_name
        self.billing_mode = billing_mode
        self.throughput = throughput
        self.key_schema = key_schema
        self.global_indexes = global_indexes
        self.attrs = attrs
        self.compression_type = compression_type

        self.failure_code: Optional[str] = None
        self.failure_message: Optional[str] = None
        self.table: Optional["Table"] = None
        self.table_arn = f"arn:{get_partition(self.region_name)}:dynamodb:{self.region_name}:{self.account_id}:table/{table_name}"

        self.processed_count = 0
        self.processed_bytes = 0
        self.error_count = 0
        self.imported_count = 0

    def run(self) -> None:
        s3_bucket_name = self.s3_source["S3Bucket"]

        try:
            from moto.s3.models import s3_backends

            s3_backend = s3_backends[self.account_id][self.partition]
            bucket = s3_backend.buckets[s3_bucket_name]
        except KeyError:
            self.status = "FAILED"
            self.failure_code = "S3NoSuchBucket"
            self.failure_message = "The specified bucket does not exist"
            return

        try:
            self._process_s3_files(s3_backend, bucket)
        except Exception as e:
            self.status = "FAILED"
            self.failure_code = "UNKNOWN"
            self.failure_message = str(e)

    def _process_s3_files(self, s3_backend: "S3Backend", bucket: "FakeBucket") -> None:
        # CREATE TABLE
        from moto.dynamodb.models import dynamodb_backends

        dynamo: DynamoDBBackend = dynamodb_backends[self.account_id][self.region_name]
        self.table = dynamo.create_table(
            name=self.table_name,
            billing_mode=self.billing_mode,
            throughput=self.throughput,
            schema=self.key_schema,
            global_indexes=self.global_indexes,
            indexes=None,
            attr=self.attrs,
            sse_specification=None,
            streams=None,
            tags=[],
            deletion_protection_enabled=False,
        )

        # Load data from S3
        keys, _, _, _ = s3_backend.list_objects(
            bucket,
            prefix=self.s3_source.get("S3KeyPrefix"),
            delimiter=None,
            marker=None,
            max_keys=None,
        )

        from py_partiql_parser import JsonParser

        for key in keys:
            if self.compression_type == "GZIP":
                content = gzip_decompress(key.value).decode("utf-8")
            else:
                content = key.value.decode("utf-8")
            result = JsonParser.parse(original=content)

            for json_object in result:
                try:
                    self.processed_count += 1
                    self.processed_bytes += len(json_object)
                    self.table.put_item(item_attrs=json_object["Item"])
                    self.imported_count += 1
                except Exception as e:
                    self.failure_message = str(e)
                    self.error_count += 1

        self.status = "COMPLETED" if self.error_count == 0 else "FAILED"

    def response(self) -> Dict[str, Any]:
        return {
            "ImportArn": self.arn,
            "ImportStatus": self.status,
            "TableArn": self.table_arn,
            "FailureCode": self.failure_code,
            "FailureMessage": self.failure_message,
            "ProcessedItemCount": self.processed_count,
            "ProcessedSizeBytes": self.processed_bytes,
            "ErrorCount": self.error_count,
            "ImportedItemCount": self.imported_count,
        }
