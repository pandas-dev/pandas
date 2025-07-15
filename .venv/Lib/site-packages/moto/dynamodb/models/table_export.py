import json
from threading import Thread
from typing import TYPE_CHECKING, Any, Dict, Optional
from uuid import uuid4

from moto.core.utils import gzip_compress
from moto.utilities.utils import get_partition

if TYPE_CHECKING:
    from moto.dynamodb.models import Table
    from moto.s3.models import S3Backend


class TableExport(Thread):
    def __init__(
        self,
        s3_bucket: str,
        s3_prefix: str,
        region_name: str,
        account_id: str,
        table: "Table",
        export_format: str,
        export_type: str,
    ):
        super().__init__()
        self.partition = get_partition(region_name)
        self.table = table
        self.arn = f"arn:{self.partition}:dynamodb:{region_name}:{account_id}:table/{table.table_arn}/import/{str(uuid4()).replace('-', '')}"
        self.s3_bucket = s3_bucket
        self.s3_prefix = s3_prefix
        self.status = "IN_PROGRESS"
        self.export_format = export_format
        self.export_type = export_type
        self.account_id = account_id
        self.region_name = region_name

        self.failure_code: Optional[str] = None
        self.failure_message: Optional[str] = None
        self.item_count = 0
        self.processed_bytes = 0
        self.error_count = 0

    def run(self) -> None:
        try:
            from moto.s3.models import s3_backends

            s3_backend = s3_backends[self.account_id][self.partition]
            s3_backend.buckets[self.s3_bucket]

        except KeyError:
            self.status = "FAILED"
            self.failure_code = "S3NoSuchBucket"
            self.failure_message = "The specified bucket does not exist"
            return

        try:
            self._backup_to_s3_file(s3_backend)
        except Exception as e:
            self.status = "FAILED"
            self.failure_code = "UNKNOWN"
            self.failure_message = str(e)

    def _backup_to_s3_file(self, s3_backend: "S3Backend") -> None:
        backup = ""
        for item in self.table.all_items():
            json_item = item.to_json(root_attr_name="Item")
            backup += json.dumps(json_item) + "\n"
            self.processed_bytes += len(json_item)
            self.item_count += 1
        content = gzip_compress(backup.encode("utf-8"))
        s3_backend.put_object(
            bucket_name=self.s3_bucket,
            key_name=f"{self.s3_prefix}/AWSDynamoDB/{str(uuid4())}/data/{str(uuid4())}.gz",
            value=content,
        )

        self.status = "COMPLETED" if self.error_count == 0 else "FAILED"

    def response(self) -> Dict[str, Any]:
        return {
            "ExportArn": self.arn,
            "ExportStatus": self.status,
            "FailureCode": self.failure_code,
            "FailureMessage": self.failure_message,
            "ExportFormat": self.export_format,
            "ExportType": self.export_type,
            "ItemCount": self.item_count,
            "BilledSizeBytes": self.processed_bytes,
        }
