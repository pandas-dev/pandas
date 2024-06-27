from collections import OrderedDict
from typing import Any, Dict, List


def cfn_to_api_encryption(
    bucket_encryption_properties: Dict[str, Any],
) -> Dict[str, Any]:
    sse_algorithm = bucket_encryption_properties["ServerSideEncryptionConfiguration"][
        0
    ]["ServerSideEncryptionByDefault"]["SSEAlgorithm"]
    kms_master_key_id = bucket_encryption_properties[
        "ServerSideEncryptionConfiguration"
    ][0]["ServerSideEncryptionByDefault"].get("KMSMasterKeyID")
    apply_server_side_encryption_by_default = OrderedDict()
    apply_server_side_encryption_by_default["SSEAlgorithm"] = sse_algorithm
    if kms_master_key_id:
        apply_server_side_encryption_by_default["KMSMasterKeyID"] = kms_master_key_id
    rule = OrderedDict(
        {"ApplyServerSideEncryptionByDefault": apply_server_side_encryption_by_default}
    )
    return OrderedDict(
        {"@xmlns": "http://s3.amazonaws.com/doc/2006-03-01/", "Rule": rule}
    )


def is_replacement_update(properties: List[str]) -> bool:
    properties_requiring_replacement_update = ["BucketName", "ObjectLockEnabled"]
    return any(
        [
            property_requiring_replacement in properties
            for property_requiring_replacement in properties_requiring_replacement_update
        ]
    )
