from moto.core.base_backend import BackendDict, BaseBackend
from moto.utilities.utils import PARTITION_NAMES


class InstanceMetadataBackend(BaseBackend):
    pass


instance_metadata_backends = BackendDict(
    InstanceMetadataBackend,
    "instance_metadata",
    use_boto3_regions=False,
    additional_regions=PARTITION_NAMES,
)
