"""EBSBackend class with methods for supported APIs."""

from typing import Any, Dict, List, Optional, Tuple

from moto.core.base_backend import BackendDict, BaseBackend
from moto.core.common_models import BaseModel
from moto.core.utils import unix_time
from moto.ec2.models import EC2Backend, ec2_backends
from moto.ec2.models.elastic_block_store import Snapshot
from moto.moto_api._internal import mock_random


class Block(BaseModel):
    def __init__(
        self, block_data: str, checksum: str, checksum_algorithm: str, data_length: str
    ):
        self.block_data = block_data
        self.checksum = checksum
        self.checksum_algorithm = checksum_algorithm
        self.data_length = data_length
        self.block_token = str(mock_random.uuid4())


class EBSSnapshot(BaseModel):
    def __init__(self, account_id: str, snapshot: Snapshot):
        self.account_id = account_id
        self.snapshot_id = snapshot.id
        self.status = "pending"
        self.start_time = unix_time()
        self.volume_size = snapshot.volume.size
        self.block_size = 512
        self.tags = [
            {"Key": t["key"], "Value": t["value"]} for t in snapshot.get_tags()
        ]
        self.description = snapshot.description

        self.blocks: Dict[str, Block] = dict()

    def put_block(
        self,
        block_idx: str,
        block_data: str,
        checksum: str,
        checksum_algorithm: str,
        data_length: str,
    ) -> None:
        block = Block(block_data, checksum, checksum_algorithm, data_length)
        self.blocks[block_idx] = block

    def to_json(self) -> Dict[str, Any]:
        return {
            "SnapshotId": self.snapshot_id,
            "OwnerId": self.account_id,
            "Status": self.status,
            "StartTime": self.start_time,
            "VolumeSize": self.volume_size,
            "BlockSize": self.block_size,
            "Tags": self.tags,
            "Description": self.description,
        }


class EBSBackend(BaseBackend):
    """Implementation of EBS APIs."""

    def __init__(self, region_name: str, account_id: str):
        super().__init__(region_name, account_id)
        self.snapshots: Dict[str, EBSSnapshot] = dict()

    @property
    def ec2_backend(self) -> EC2Backend:
        return ec2_backends[self.account_id][self.region_name]

    def start_snapshot(
        self, volume_size: int, tags: Optional[List[Dict[str, str]]], description: str
    ) -> EBSSnapshot:
        zone_name = f"{self.region_name}a"
        vol = self.ec2_backend.create_volume(size=volume_size, zone_name=zone_name)
        snapshot = self.ec2_backend.create_snapshot(
            volume_id=vol.id, description=description
        )
        if tags:
            snapshot.add_tags({tag["Key"]: tag["Value"] for tag in tags})
        ebs_snapshot = EBSSnapshot(account_id=self.account_id, snapshot=snapshot)
        self.snapshots[ebs_snapshot.snapshot_id] = ebs_snapshot
        return ebs_snapshot

    def complete_snapshot(self, snapshot_id: str) -> Dict[str, str]:
        self.snapshots[snapshot_id].status = "completed"
        return {"Status": "completed"}

    def put_snapshot_block(
        self,
        snapshot_id: str,
        block_index: str,
        block_data: str,
        checksum: str,
        checksum_algorithm: str,
        data_length: str,
    ) -> Tuple[str, str]:
        snapshot = self.snapshots[snapshot_id]
        snapshot.put_block(
            block_index, block_data, checksum, checksum_algorithm, data_length
        )
        return checksum, checksum_algorithm

    def get_snapshot_block(self, snapshot_id: str, block_index: str) -> Block:
        """
        The BlockToken-parameter is not yet implemented
        """
        snapshot = self.snapshots[snapshot_id]
        return snapshot.blocks[block_index]

    def list_changed_blocks(
        self, first_snapshot_id: str, second_snapshot_id: str
    ) -> Tuple[Dict[str, Tuple[str, Optional[str]]], EBSSnapshot]:
        """
        The following parameters are not yet implemented: NextToken, MaxResults, StartingBlockIndex
        """
        snapshot1 = self.snapshots[first_snapshot_id]
        snapshot2 = self.snapshots[second_snapshot_id]
        changed_blocks: Dict[str, Tuple[str, Optional[str]]] = (
            dict()
        )  # {idx: (token1, token2), ..}
        for idx in snapshot1.blocks:
            block1 = snapshot1.blocks[idx]
            if idx in snapshot2.blocks:
                block2 = snapshot2.blocks[idx]
                if block1.block_data != block2.block_data:
                    changed_blocks[idx] = (block1.block_token, block2.block_token)
            else:
                changed_blocks[idx] = (block1.block_token, None)

        return changed_blocks, snapshot1

    def list_snapshot_blocks(self, snapshot_id: str) -> EBSSnapshot:
        return self.snapshots[snapshot_id]


ebs_backends = BackendDict(EBSBackend, "ebs")
