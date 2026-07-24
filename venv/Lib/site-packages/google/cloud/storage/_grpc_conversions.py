# Copyright 2026 Google LLC
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

from google.protobuf import timestamp_pb2

from google.cloud import _storage_v2

# Map Python Blob attributes to GCS V2 Object proto field names.
_BLOB_ATTR_TO_PROTO_FIELD = {
    "content_type": "content_type",
    "metadata": "metadata",
    "kms_key_name": "kms_key",
    "cache_control": "cache_control",
    "content_disposition": "content_disposition",
    "content_encoding": "content_encoding",
    "content_language": "content_language",
    "temporary_hold": "temporary_hold",
    "event_based_hold": "event_based_hold",
}


def blob_to_proto(blob):
    """Converts a Blob instance to a GCS V2 Object proto message."""

    resource_params = {
        "name": blob.name,
    }

    if blob.bucket:
        resource_params["bucket"] = f"projects/_/buckets/{blob.bucket.name}"

    for attr_name, proto_field in _BLOB_ATTR_TO_PROTO_FIELD.items():
        value = getattr(blob, attr_name, None)
        if value is not None:
            resource_params[proto_field] = value

    custom_time = getattr(blob, "custom_time", None)
    if custom_time is not None:
        custom_time_proto = timestamp_pb2.Timestamp()
        custom_time_proto.FromDatetime(custom_time)
        resource_params["custom_time"] = custom_time_proto

    acl = getattr(blob, "acl", None)
    if acl is not None and getattr(acl, "loaded", False):
        acl_entries = []
        for entry in acl:
            acl_entries.append(
                _storage_v2.ObjectAccessControl(
                    role=entry["role"],
                    entity=entry["entity"],
                )
            )
        if acl_entries:
            resource_params["acl"] = acl_entries

    retention = getattr(blob, "retention", None)
    if retention:
        mode_str = retention.get("mode")
        mode = _storage_v2.Object.Retention.Mode.MODE_UNSPECIFIED
        if mode_str:
            # GCS retention modes are 'Locked' or 'Unlocked'
            mode = getattr(
                _storage_v2.Object.Retention.Mode,
                mode_str.upper(),
                _storage_v2.Object.Retention.Mode.MODE_UNSPECIFIED,
            )

        retain_until_time_proto = None
        retain_until_time = retention.get("retain_until_time")
        if retain_until_time is not None:
            retain_until_time_proto = timestamp_pb2.Timestamp()
            retain_until_time_proto.FromDatetime(retain_until_time)

        resource_params["retention"] = _storage_v2.Object.Retention(
            mode=mode,
            retain_until_time=retain_until_time_proto,
        )

    contexts = getattr(blob, "contexts", None)
    if contexts:
        custom_contexts = {}
        for key, payload in contexts.custom.items():
            custom_contexts[key] = _storage_v2.ObjectCustomContextPayload(
                value=payload.value
            )

        resource_params["contexts"] = _storage_v2.ObjectContexts(custom=custom_contexts)

    return _storage_v2.Object(**resource_params)
