# Copyright 2025 Google LLC
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

from typing import ClassVar as _ClassVar
from typing import Iterable as _Iterable
from typing import Mapping as _Mapping
from typing import Optional as _Optional
from typing import Union as _Union

from google.api import annotations_pb2 as _annotations_pb2
from google.api import client_pb2 as _client_pb2
from google.api import field_behavior_pb2 as _field_behavior_pb2
from google.api import resource_pb2 as _resource_pb2
from google.iam.v1 import options_pb2 as _options_pb2
from google.iam.v1 import policy_pb2 as _policy_pb2
from google.protobuf import descriptor as _descriptor
from google.protobuf import field_mask_pb2 as _field_mask_pb2
from google.protobuf import message as _message
from google.protobuf.internal import containers as _containers

DESCRIPTOR: _descriptor.FileDescriptor

class SetIamPolicyRequest(_message.Message):
    __slots__ = ("resource", "policy", "update_mask")
    RESOURCE_FIELD_NUMBER: _ClassVar[int]
    POLICY_FIELD_NUMBER: _ClassVar[int]
    UPDATE_MASK_FIELD_NUMBER: _ClassVar[int]
    resource: str
    policy: _policy_pb2.Policy
    update_mask: _field_mask_pb2.FieldMask
    def __init__(
        self,
        resource: _Optional[str] = ...,
        policy: _Optional[_Union[_policy_pb2.Policy, _Mapping]] = ...,
        update_mask: _Optional[_Union[_field_mask_pb2.FieldMask, _Mapping]] = ...,
    ) -> None: ...

class GetIamPolicyRequest(_message.Message):
    __slots__ = ("resource", "options")
    RESOURCE_FIELD_NUMBER: _ClassVar[int]
    OPTIONS_FIELD_NUMBER: _ClassVar[int]
    resource: str
    options: _options_pb2.GetPolicyOptions
    def __init__(
        self,
        resource: _Optional[str] = ...,
        options: _Optional[_Union[_options_pb2.GetPolicyOptions, _Mapping]] = ...,
    ) -> None: ...

class TestIamPermissionsRequest(_message.Message):
    __slots__ = ("resource", "permissions")
    RESOURCE_FIELD_NUMBER: _ClassVar[int]
    PERMISSIONS_FIELD_NUMBER: _ClassVar[int]
    resource: str
    permissions: _containers.RepeatedScalarFieldContainer[str]
    def __init__(
        self,
        resource: _Optional[str] = ...,
        permissions: _Optional[_Iterable[str]] = ...,
    ) -> None: ...

class TestIamPermissionsResponse(_message.Message):
    __slots__ = ("permissions",)
    PERMISSIONS_FIELD_NUMBER: _ClassVar[int]
    permissions: _containers.RepeatedScalarFieldContainer[str]
    def __init__(self, permissions: _Optional[_Iterable[str]] = ...) -> None: ...
