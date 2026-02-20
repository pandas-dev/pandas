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

from google.protobuf import descriptor as _descriptor
from google.protobuf import message as _message
from google.protobuf.internal import containers as _containers
from google.protobuf.internal import enum_type_wrapper as _enum_type_wrapper
from google.type import expr_pb2 as _expr_pb2

DESCRIPTOR: _descriptor.FileDescriptor

class Policy(_message.Message):
    __slots__ = ("version", "bindings", "audit_configs", "etag")
    VERSION_FIELD_NUMBER: _ClassVar[int]
    BINDINGS_FIELD_NUMBER: _ClassVar[int]
    AUDIT_CONFIGS_FIELD_NUMBER: _ClassVar[int]
    ETAG_FIELD_NUMBER: _ClassVar[int]
    version: int
    bindings: _containers.RepeatedCompositeFieldContainer[Binding]
    audit_configs: _containers.RepeatedCompositeFieldContainer[AuditConfig]
    etag: bytes
    def __init__(
        self,
        version: _Optional[int] = ...,
        bindings: _Optional[_Iterable[_Union[Binding, _Mapping]]] = ...,
        audit_configs: _Optional[_Iterable[_Union[AuditConfig, _Mapping]]] = ...,
        etag: _Optional[bytes] = ...,
    ) -> None: ...

class Binding(_message.Message):
    __slots__ = ("role", "members", "condition")
    ROLE_FIELD_NUMBER: _ClassVar[int]
    MEMBERS_FIELD_NUMBER: _ClassVar[int]
    CONDITION_FIELD_NUMBER: _ClassVar[int]
    role: str
    members: _containers.RepeatedScalarFieldContainer[str]
    condition: _expr_pb2.Expr
    def __init__(
        self,
        role: _Optional[str] = ...,
        members: _Optional[_Iterable[str]] = ...,
        condition: _Optional[_Union[_expr_pb2.Expr, _Mapping]] = ...,
    ) -> None: ...

class AuditConfig(_message.Message):
    __slots__ = ("service", "audit_log_configs")
    SERVICE_FIELD_NUMBER: _ClassVar[int]
    AUDIT_LOG_CONFIGS_FIELD_NUMBER: _ClassVar[int]
    service: str
    audit_log_configs: _containers.RepeatedCompositeFieldContainer[AuditLogConfig]
    def __init__(
        self,
        service: _Optional[str] = ...,
        audit_log_configs: _Optional[_Iterable[_Union[AuditLogConfig, _Mapping]]] = ...,
    ) -> None: ...

class AuditLogConfig(_message.Message):
    __slots__ = ("log_type", "exempted_members")

    class LogType(int, metaclass=_enum_type_wrapper.EnumTypeWrapper):
        __slots__ = ()
        LOG_TYPE_UNSPECIFIED: _ClassVar[AuditLogConfig.LogType]
        ADMIN_READ: _ClassVar[AuditLogConfig.LogType]
        DATA_WRITE: _ClassVar[AuditLogConfig.LogType]
        DATA_READ: _ClassVar[AuditLogConfig.LogType]
    LOG_TYPE_UNSPECIFIED: AuditLogConfig.LogType
    ADMIN_READ: AuditLogConfig.LogType
    DATA_WRITE: AuditLogConfig.LogType
    DATA_READ: AuditLogConfig.LogType
    LOG_TYPE_FIELD_NUMBER: _ClassVar[int]
    EXEMPTED_MEMBERS_FIELD_NUMBER: _ClassVar[int]
    log_type: AuditLogConfig.LogType
    exempted_members: _containers.RepeatedScalarFieldContainer[str]
    def __init__(
        self,
        log_type: _Optional[_Union[AuditLogConfig.LogType, str]] = ...,
        exempted_members: _Optional[_Iterable[str]] = ...,
    ) -> None: ...

class PolicyDelta(_message.Message):
    __slots__ = ("binding_deltas", "audit_config_deltas")
    BINDING_DELTAS_FIELD_NUMBER: _ClassVar[int]
    AUDIT_CONFIG_DELTAS_FIELD_NUMBER: _ClassVar[int]
    binding_deltas: _containers.RepeatedCompositeFieldContainer[BindingDelta]
    audit_config_deltas: _containers.RepeatedCompositeFieldContainer[AuditConfigDelta]
    def __init__(
        self,
        binding_deltas: _Optional[_Iterable[_Union[BindingDelta, _Mapping]]] = ...,
        audit_config_deltas: _Optional[
            _Iterable[_Union[AuditConfigDelta, _Mapping]]
        ] = ...,
    ) -> None: ...

class BindingDelta(_message.Message):
    __slots__ = ("action", "role", "member", "condition")

    class Action(int, metaclass=_enum_type_wrapper.EnumTypeWrapper):
        __slots__ = ()
        ACTION_UNSPECIFIED: _ClassVar[BindingDelta.Action]
        ADD: _ClassVar[BindingDelta.Action]
        REMOVE: _ClassVar[BindingDelta.Action]
    ACTION_UNSPECIFIED: BindingDelta.Action
    ADD: BindingDelta.Action
    REMOVE: BindingDelta.Action
    ACTION_FIELD_NUMBER: _ClassVar[int]
    ROLE_FIELD_NUMBER: _ClassVar[int]
    MEMBER_FIELD_NUMBER: _ClassVar[int]
    CONDITION_FIELD_NUMBER: _ClassVar[int]
    action: BindingDelta.Action
    role: str
    member: str
    condition: _expr_pb2.Expr
    def __init__(
        self,
        action: _Optional[_Union[BindingDelta.Action, str]] = ...,
        role: _Optional[str] = ...,
        member: _Optional[str] = ...,
        condition: _Optional[_Union[_expr_pb2.Expr, _Mapping]] = ...,
    ) -> None: ...

class AuditConfigDelta(_message.Message):
    __slots__ = ("action", "service", "exempted_member", "log_type")

    class Action(int, metaclass=_enum_type_wrapper.EnumTypeWrapper):
        __slots__ = ()
        ACTION_UNSPECIFIED: _ClassVar[AuditConfigDelta.Action]
        ADD: _ClassVar[AuditConfigDelta.Action]
        REMOVE: _ClassVar[AuditConfigDelta.Action]
    ACTION_UNSPECIFIED: AuditConfigDelta.Action
    ADD: AuditConfigDelta.Action
    REMOVE: AuditConfigDelta.Action
    ACTION_FIELD_NUMBER: _ClassVar[int]
    SERVICE_FIELD_NUMBER: _ClassVar[int]
    EXEMPTED_MEMBER_FIELD_NUMBER: _ClassVar[int]
    LOG_TYPE_FIELD_NUMBER: _ClassVar[int]
    action: AuditConfigDelta.Action
    service: str
    exempted_member: str
    log_type: str
    def __init__(
        self,
        action: _Optional[_Union[AuditConfigDelta.Action, str]] = ...,
        service: _Optional[str] = ...,
        exempted_member: _Optional[str] = ...,
        log_type: _Optional[str] = ...,
    ) -> None: ...
