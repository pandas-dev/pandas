"""Handles incoming ebs requests, invokes methods, returns responses."""

import json
from typing import Any

from moto.core.common_types import TYPE_RESPONSE
from moto.core.responses import BaseResponse

from .models import EBSBackend, ebs_backends


class EBSResponse(BaseResponse):
    """Handler for EBS requests and responses."""

    def __init__(self) -> None:
        super().__init__(service_name="ebs")

    def setup_class(self, request: Any, full_url: str, headers: Any) -> None:  # type: ignore
        super().setup_class(request, full_url, headers, use_raw_body=True)

    @property
    def ebs_backend(self) -> EBSBackend:
        """Return backend instance specific for this region."""
        return ebs_backends[self.current_account][self.region]

    def start_snapshot(self) -> str:
        params = json.loads(self.body)
        volume_size = params.get("VolumeSize")
        tags = params.get("Tags")
        description = params.get("Description")
        snapshot = self.ebs_backend.start_snapshot(
            volume_size=volume_size,
            tags=tags,
            description=description,
        )
        return json.dumps(snapshot.to_json())

    def complete_snapshot(self) -> TYPE_RESPONSE:
        snapshot_id = self.parsed_url.path.split("/")[-1]
        status = self.ebs_backend.complete_snapshot(snapshot_id=snapshot_id)
        return 202, {"status": 202}, json.dumps(status)

    def put_snapshot_block(self) -> TYPE_RESPONSE:
        snapshot_id = self.parsed_url.path.split("/")[-3]
        block_index = self.parsed_url.path.split("/")[-1]
        block_data = self.body
        headers = {k.lower(): v for k, v in self.headers.items()}
        checksum = headers.get("x-amz-checksum")
        checksum_algorithm = headers.get("x-amz-checksum-algorithm")
        data_length = headers.get("x-amz-data-length")
        checksum, checksum_algorithm = self.ebs_backend.put_snapshot_block(
            snapshot_id=snapshot_id,
            block_index=block_index,
            block_data=block_data,
            checksum=checksum,  # type: ignore
            checksum_algorithm=checksum_algorithm,  # type: ignore
            data_length=data_length,  # type: ignore
        )
        resp_headers = {
            "status": 201,
            "x-amz-Checksum": checksum,
            "x-amz-Checksum-Algorithm": checksum_algorithm,
        }
        return 201, resp_headers, "{}"

    def get_snapshot_block(self) -> TYPE_RESPONSE:
        snapshot_id = self.path.split("/")[-3]
        block_index = self.path.split("/")[-1]
        block = self.ebs_backend.get_snapshot_block(
            snapshot_id=snapshot_id,
            block_index=block_index,
        )
        headers = {
            "x-amz-Checksum": block.checksum,
            "x-amz-Checksum-Algorithm": block.checksum_algorithm,
            "x-amz-Data-Length": block.data_length,
        }
        return 200, headers, block.block_data

    def list_changed_blocks(self) -> str:
        first_snapshot_id = self._get_params().get("firstSnapshotId")
        second_snapshot_id = self.path.split("/")[-2]
        changed_blocks, snapshot = self.ebs_backend.list_changed_blocks(
            first_snapshot_id=first_snapshot_id,  # type: ignore[arg-type]
            second_snapshot_id=second_snapshot_id,
        )
        blocks = [
            {"BlockIndex": idx, "FirstBlockToken": x, "SecondBlockToken": y}
            for idx, (x, y) in changed_blocks.items()
        ]
        return json.dumps(
            dict(
                ChangedBlocks=blocks,
                VolumeSize=snapshot.volume_size,
                BlockSize=snapshot.block_size,
            )
        )

    def list_snapshot_blocks(self) -> str:
        snapshot_id = self.path.split("/")[-2]
        snapshot = self.ebs_backend.list_snapshot_blocks(
            snapshot_id=snapshot_id,
        )
        blocks = [
            {"BlockIndex": idx, "BlockToken": b.block_token}
            for idx, b in snapshot.blocks.items()
        ]
        return json.dumps(
            {
                "Blocks": blocks,
                "VolumeSize": snapshot.volume_size,
                "BlockSize": snapshot.block_size,
            }
        )
