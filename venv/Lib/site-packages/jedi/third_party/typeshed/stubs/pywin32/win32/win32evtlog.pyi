from _typeshed import Incomplete
from collections.abc import Iterable

import _win32typing
from win32.lib.pywintypes import error as error

def ReadEventLog(
    Handle: _win32typing.PyEVTLOG_HANDLE, Flags: int, Offset: int, Size=..., /
) -> list[_win32typing.PyEventLogRecord]: ...
def ClearEventLog(handle: _win32typing.PyEVTLOG_HANDLE, eventLogName: str, /) -> None: ...
def BackupEventLog(handle, eventLogName: str, /) -> None: ...
def CloseEventLog(handle: _win32typing.PyEVTLOG_HANDLE, /) -> None: ...
def DeregisterEventSource(handle, /) -> None: ...
def NotifyChangeEventLog(handle, handle1, /) -> None: ...
def GetNumberOfEventLogRecords(handle: _win32typing.PyEVTLOG_HANDLE, /) -> int: ...
def GetOldestEventLogRecord(): ...
def OpenEventLog(serverName: str | None, sourceName: str, /) -> _win32typing.PyEVTLOG_HANDLE: ...
def RegisterEventSource(serverName: str | None, sourceName: str, /): ...
def OpenBackupEventLog(serverName: str, fileName: str, /) -> _win32typing.PyEVTLOG_HANDLE: ...
def ReportEvent(
    EventLog: int,
    Type: int,
    Category: int,
    EventID: int,
    UserSid: _win32typing.PySID | None,
    Strings: Iterable[str] | None,
    RawData: bytes | None,
    /,
) -> None: ...
def EvtOpenChannelEnum(Session: _win32typing.PyEVT_HANDLE | None = ..., Flags: int = ...) -> _win32typing.PyEVT_HANDLE: ...
def EvtFormatMessage(Metadata, Event, Flags, ResourceId=...): ...
def EvtNextChannelPath(ChannelEnum: _win32typing.PyEVT_HANDLE): ...
def EvtOpenLog(Path, Flags, Session: _win32typing.PyEVT_HANDLE | None = ...) -> _win32typing.PyEVT_HANDLE: ...
def EvtClearLog(
    ChannelPath, TargetFilePath: Incomplete | None = ..., Session: _win32typing.PyEVT_HANDLE | None = ..., Flags: int = ...
) -> None: ...
def EvtExportLog(
    Path, TargetFilePath, Flags, Query: Incomplete | None = ..., Session: _win32typing.PyEVT_HANDLE | None = ...
) -> None: ...
def EvtArchiveExportedLog(LogFilePath, Locale, Session: _win32typing.PyEVT_HANDLE | None = ..., Flags=...) -> None: ...
def EvtGetExtendedStatus(): ...
def EvtQuery(
    Path: str, Flags: int, Query: str | None = ..., Session: _win32typing.PyEVT_HANDLE | None = ...
) -> _win32typing.PyEVT_HANDLE: ...
def EvtNext(
    ResultSet: _win32typing.PyEVT_HANDLE, Count: int, Timeout: int = ..., Flags: int = ...
) -> tuple[_win32typing.PyEVT_HANDLE, ...]: ...
def EvtSeek(
    ResultSet: _win32typing.PyEVT_HANDLE, Position, Flags, Bookmark: _win32typing.PyEVT_HANDLE | None = ..., Timeout: int = ...
) -> None: ...
def EvtCreateRenderContext(Flags): ...
def EvtRender(Event: _win32typing.PyEVT_HANDLE, Flags: int, Context=...): ...
def EvtSubscribe(
    ChannelPath,
    Flags,
    SignalEvent: Incomplete | None = ...,
    Callback: Incomplete | None = ...,
    Context: Incomplete | None = ...,
    Query: Incomplete | None = ...,
    Session: _win32typing.PyEVT_HANDLE | None = ...,
    Bookmark: _win32typing.PyEVT_HANDLE | None = ...,
) -> _win32typing.PyEVT_HANDLE: ...
def EvtCreateBookmark(BookmarkXML: Incomplete | None = ...) -> _win32typing.PyEVT_HANDLE: ...
def EvtUpdateBookmark(Bookmark: _win32typing.PyEVT_HANDLE, Event: _win32typing.PyEVT_HANDLE) -> _win32typing.PyEVT_HANDLE: ...
def EvtGetChannelConfigProperty(
    ChannelConfig: _win32typing.PyEVT_HANDLE, PropertyId, Flags=...
) -> tuple[Incomplete, Incomplete]: ...
def EvtOpenChannelConfig(
    ChannelPath, Session: _win32typing.PyEVT_HANDLE | None = ..., Flags=...
) -> _win32typing.PyEVT_HANDLE: ...
def EvtOpenSession(
    Login: _win32typing.PyEVT_RPC_LOGIN, LoginClass, Timeout: int = ..., Flags=...
) -> _win32typing.PyEVT_HANDLE: ...
def EvtOpenPublisherEnum(Session: _win32typing.PyEVT_HANDLE | None = ..., Flags: int = ...) -> _win32typing.PyEVT_HANDLE: ...
def EvtNextPublisherId(PublisherEnum: _win32typing.PyEVT_HANDLE): ...
def EvtOpenPublisherMetadata(
    PublisherIdentity: str,
    Session: _win32typing.PyEVT_HANDLE | None = ...,
    LogFilePath: Incomplete | None = ...,
    Locale: int = ...,
    Flags: int = ...,
) -> _win32typing.PyEVT_HANDLE: ...
def EvtGetPublisherMetadataProperty(
    PublisherMetadata: _win32typing.PyEVT_HANDLE, PropertyId, Flags=...
) -> tuple[Incomplete, Incomplete]: ...
def EvtOpenEventMetadataEnum(PublisherMetadata: _win32typing.PyEVT_HANDLE, Flags=...) -> _win32typing.PyEVT_HANDLE: ...
def EvtNextEventMetadata(EventMetadataEnum: _win32typing.PyEVT_HANDLE, Flags=...) -> _win32typing.PyEVT_HANDLE: ...
def EvtGetEventMetadataProperty(
    EventMetadata: _win32typing.PyEVT_HANDLE, PropertyId, Flags=...
) -> tuple[Incomplete, Incomplete]: ...
def EvtGetLogInfo(Log: _win32typing.PyEVT_HANDLE, PropertyId) -> tuple[Incomplete, Incomplete]: ...
def EvtGetEventInfo(Event: _win32typing.PyEVT_HANDLE, PropertyId) -> tuple[Incomplete, Incomplete]: ...
def EvtGetObjectArraySize(ObjectArray: _win32typing.PyEVT_HANDLE): ...
def EvtGetObjectArrayProperty(
    ObjectArray: _win32typing.PyEVT_HANDLE, PropertyId, ArrayIndex, Flags=...
) -> tuple[Incomplete, Incomplete]: ...

EVENTLOG_AUDIT_FAILURE: int
EVENTLOG_AUDIT_SUCCESS: int
EVENTLOG_BACKWARDS_READ: int
EVENTLOG_END_ALL_PAIRED_EVENTS: int
EVENTLOG_END_PAIRED_EVENT: int
EVENTLOG_ERROR_TYPE: int
EVENTLOG_FORWARDS_READ: int
EVENTLOG_INFORMATION_TYPE: int
EVENTLOG_PAIRED_EVENT_ACTIVE: int
EVENTLOG_PAIRED_EVENT_INACTIVE: int
EVENTLOG_SEEK_READ: int
EVENTLOG_SEQUENTIAL_READ: int
EVENTLOG_START_PAIRED_EVENT: int
EVENTLOG_SUCCESS: int
EVENTLOG_WARNING_TYPE: int
EventMetadataEventChannel: int
EventMetadataEventID: int
EventMetadataEventKeyword: int
EventMetadataEventLevel: int
EventMetadataEventMessageID: int
EventMetadataEventOpcode: int
EventMetadataEventTask: int
EventMetadataEventTemplate: int
EventMetadataEventVersion: int
EvtChannelConfigAccess: int
EvtChannelConfigClassicEventlog: int
EvtChannelConfigEnabled: int
EvtChannelConfigIsolation: int
EvtChannelConfigOwningPublisher: int
EvtChannelConfigPropertyIdEND: int
EvtChannelConfigType: int
EvtChannelLoggingConfigAutoBackup: int
EvtChannelLoggingConfigLogFilePath: int
EvtChannelLoggingConfigMaxSize: int
EvtChannelLoggingConfigRetention: int
EvtChannelPublishingConfigBufferSize: int
EvtChannelPublishingConfigClockType: int
EvtChannelPublishingConfigControlGuid: int
EvtChannelPublishingConfigKeywords: int
EvtChannelPublishingConfigLatency: int
EvtChannelPublishingConfigLevel: int
EvtChannelPublishingConfigMaxBuffers: int
EvtChannelPublishingConfigMinBuffers: int
EvtChannelPublishingConfigSidType: int
EvtEventMetadataPropertyIdEND: int
EvtEventPath: int
EvtEventPropertyIdEND: int
EvtEventQueryIDs: int
EvtExportLogChannelPath: int
EvtExportLogFilePath: int
EvtExportLogTolerateQueryErrors: int
EvtLogAttributes: int
EvtLogCreationTime: int
EvtLogFileSize: int
EvtLogFull: int
EvtLogLastAccessTime: int
EvtLogLastWriteTime: int
EvtLogNumberOfLogRecords: int
EvtLogOldestRecordNumber: int
EvtOpenChannelPath: int
EvtOpenFilePath: int
EvtPublisherMetadataChannelReferenceFlags: int
EvtPublisherMetadataChannelReferenceID: int
EvtPublisherMetadataChannelReferenceIndex: int
EvtPublisherMetadataChannelReferenceMessageID: int
EvtPublisherMetadataChannelReferencePath: int
EvtPublisherMetadataChannelReferences: int
EvtPublisherMetadataHelpLink: int
EvtPublisherMetadataKeywordMessageID: int
EvtPublisherMetadataKeywordName: int
EvtPublisherMetadataKeywords: int
EvtPublisherMetadataKeywordValue: int
EvtPublisherMetadataLevelMessageID: int
EvtPublisherMetadataLevelName: int
EvtPublisherMetadataLevels: int
EvtPublisherMetadataLevelValue: int
EvtPublisherMetadataMessageFilePath: int
EvtPublisherMetadataOpcodeMessageID: int
EvtPublisherMetadataOpcodeName: int
EvtPublisherMetadataOpcodes: int
EvtPublisherMetadataOpcodeValue: int
EvtPublisherMetadataParameterFilePath: int
EvtPublisherMetadataPropertyIdEND: int
EvtPublisherMetadataPublisherGuid: int
EvtPublisherMetadataPublisherMessageID: int
EvtPublisherMetadataResourceFilePath: int
EvtPublisherMetadataTaskEventGuid: int
EvtPublisherMetadataTaskMessageID: int
EvtPublisherMetadataTaskName: int
EvtPublisherMetadataTasks: int
EvtPublisherMetadataTaskValue: int
EvtQueryChannelPath: int
EvtQueryFilePath: int
EvtQueryForwardDirection: int
EvtQueryReverseDirection: int
EvtQueryTolerateQueryErrors: int
EvtRenderBookmark: int
EvtRenderEventValues: int
EvtRenderEventXml: int
EvtRpcLogin: int
EvtRpcLoginAuthDefault: int
EvtRpcLoginAuthKerberos: int
EvtRpcLoginAuthNegotiate: int
EvtRpcLoginAuthNTLM: int
EvtSeekOriginMask: int
EvtSeekRelativeToBookmark: int
EvtSeekRelativeToCurrent: int
EvtSeekRelativeToFirst: int
EvtSeekRelativeToLast: int
EvtSeekStrict: int
EvtSubscribeActionDeliver: int
EvtSubscribeActionError: int
EvtSubscribeOriginMask: int
EvtSubscribeStartAfterBookmark: int
EvtSubscribeStartAtOldestRecord: int
EvtSubscribeStrict: int
EvtSubscribeToFutureEvents: int
EvtSubscribeTolerateQueryErrors: int
EvtVarTypeAnsiString: int
EvtVarTypeBinary: int
EvtVarTypeBoolean: int
EvtVarTypeByte: int
EvtVarTypeDouble: int
EvtVarTypeEvtHandle: int
EvtVarTypeEvtXml: int
EvtVarTypeFileTime: int
EvtVarTypeGuid: int
EvtVarTypeHexInt32: int
EvtVarTypeHexInt64: int
EvtVarTypeInt16: int
EvtVarTypeInt32: int
EvtVarTypeInt64: int
EvtVarTypeNull: int
EvtVarTypeSByte: int
EvtVarTypeSid: int
EvtVarTypeSingle: int
EvtVarTypeSizeT: int
EvtVarTypeString: int
EvtVarTypeSysTime: int
EvtVarTypeUInt16: int
EvtVarTypeUInt32: int
EvtVarTypeUInt64: int
EvtChannelPublisherList: int
EvtFormatMessageChannel: int
EvtFormatMessageEvent: int
EvtFormatMessageId: int
EvtFormatMessageKeyword: int
EvtFormatMessageLevel: int
EvtFormatMessageOpcode: int
EvtFormatMessageProvider: int
EvtFormatMessageTask: int
EvtFormatMessageXml: int
EvtRenderContextSystem: int
EvtRenderContextUser: int
EvtRenderContextValues: int
EvtSystemActivityID: int
EvtSystemChannel: int
EvtSystemComputer: int
EvtSystemEventID: int
EvtSystemEventRecordId: int
EvtSystemKeywords: int
EvtSystemLevel: int
EvtSystemOpcode: int
EvtSystemProcessID: int
EvtSystemPropertyIdEND: int
EvtSystemProviderGuid: int
EvtSystemProviderName: int
EvtSystemQualifiers: int
EvtSystemRelatedActivityID: int
EvtSystemTask: int
EvtSystemThreadID: int
EvtSystemTimeCreated: int
EvtSystemUserID: int
EvtSystemVersion: int
UNICODE: int
