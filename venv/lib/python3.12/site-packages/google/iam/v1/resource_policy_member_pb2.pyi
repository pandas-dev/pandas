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
from typing import Optional as _Optional

from google.api import field_behavior_pb2 as _field_behavior_pb2
from google.protobuf import descriptor as _descriptor
from google.protobuf import message as _message

DESCRIPTOR: _descriptor.FileDescriptor

class ResourcePolicyMember(_message.Message):
    __slots__ = ("iam_policy_name_principal", "iam_policy_uid_principal")
    IAM_POLICY_NAME_PRINCIPAL_FIELD_NUMBER: _ClassVar[int]
    IAM_POLICY_UID_PRINCIPAL_FIELD_NUMBER: _ClassVar[int]
    iam_policy_name_principal: str
    iam_policy_uid_principal: str
    def __init__(
        self,
        iam_policy_name_principal: _Optional[str] = ...,
        iam_policy_uid_principal: _Optional[str] = ...,
    ) -> None: ...
