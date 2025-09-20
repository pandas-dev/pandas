from datetime import datetime
from enum import Enum
from typing import IO, Dict, Iterable, List, Optional, TypedDict, Union

from ..exceptions import AWSError as ServiceException

AliasDescription = str
Arn = str
CharacterRestrictedName = str
ClientToken = str
ConnectorParameters = str
Definition = str
Enabled = bool
ErrorMessage = str
EvaluationFailureLocation = str
HTTPBody = str
HTTPHeaders = str
HTTPMethod = str
HTTPProtocol = str
HTTPStatusCode = str
HTTPStatusMessage = str
Identity = str
IncludeExecutionData = bool
IncludeExecutionDataGetExecutionHistory = bool
KmsDataKeyReusePeriodSeconds = int
KmsKeyId = str
ListExecutionsPageToken = str
LongArn = str
MapRunLabel = str
MaxConcurrency = int
Name = str
PageSize = int
PageToken = str
Publish = bool
RedriveCount = int
RevealSecrets = bool
ReverseOrder = bool
RevisionId = str
SensitiveCause = str
SensitiveData = str
SensitiveDataJobInput = str
SensitiveError = str
StateName = str
TagKey = str
TagValue = str
TaskToken = str
ToleratedFailurePercentage = float
TraceHeader = str
URL = str
UnsignedInteger = int
ValidateStateMachineDefinitionCode = str
ValidateStateMachineDefinitionLocation = str
ValidateStateMachineDefinitionMaxResult = int
ValidateStateMachineDefinitionMessage = str
ValidateStateMachineDefinitionTruncated = bool
VariableName = str
VariableValue = str
VersionDescription = str
VersionWeight = int
includedDetails = bool
truncated = bool


class EncryptionType(str, Enum):
    AWS_OWNED_KEY = "AWS_OWNED_KEY"
    CUSTOMER_MANAGED_KMS_KEY = "CUSTOMER_MANAGED_KMS_KEY"


class ExecutionRedriveFilter(str, Enum):
    REDRIVEN = "REDRIVEN"
    NOT_REDRIVEN = "NOT_REDRIVEN"


class ExecutionRedriveStatus(str, Enum):
    REDRIVABLE = "REDRIVABLE"
    NOT_REDRIVABLE = "NOT_REDRIVABLE"
    REDRIVABLE_BY_MAP_RUN = "REDRIVABLE_BY_MAP_RUN"


class ExecutionStatus(str, Enum):
    RUNNING = "RUNNING"
    SUCCEEDED = "SUCCEEDED"
    FAILED = "FAILED"
    TIMED_OUT = "TIMED_OUT"
    ABORTED = "ABORTED"
    PENDING_REDRIVE = "PENDING_REDRIVE"


class HistoryEventType(str, Enum):
    ActivityFailed = "ActivityFailed"
    ActivityScheduled = "ActivityScheduled"
    ActivityScheduleFailed = "ActivityScheduleFailed"
    ActivityStarted = "ActivityStarted"
    ActivitySucceeded = "ActivitySucceeded"
    ActivityTimedOut = "ActivityTimedOut"
    ChoiceStateEntered = "ChoiceStateEntered"
    ChoiceStateExited = "ChoiceStateExited"
    ExecutionAborted = "ExecutionAborted"
    ExecutionFailed = "ExecutionFailed"
    ExecutionStarted = "ExecutionStarted"
    ExecutionSucceeded = "ExecutionSucceeded"
    ExecutionTimedOut = "ExecutionTimedOut"
    FailStateEntered = "FailStateEntered"
    LambdaFunctionFailed = "LambdaFunctionFailed"
    LambdaFunctionScheduled = "LambdaFunctionScheduled"
    LambdaFunctionScheduleFailed = "LambdaFunctionScheduleFailed"
    LambdaFunctionStarted = "LambdaFunctionStarted"
    LambdaFunctionStartFailed = "LambdaFunctionStartFailed"
    LambdaFunctionSucceeded = "LambdaFunctionSucceeded"
    LambdaFunctionTimedOut = "LambdaFunctionTimedOut"
    MapIterationAborted = "MapIterationAborted"
    MapIterationFailed = "MapIterationFailed"
    MapIterationStarted = "MapIterationStarted"
    MapIterationSucceeded = "MapIterationSucceeded"
    MapStateAborted = "MapStateAborted"
    MapStateEntered = "MapStateEntered"
    MapStateExited = "MapStateExited"
    MapStateFailed = "MapStateFailed"
    MapStateStarted = "MapStateStarted"
    MapStateSucceeded = "MapStateSucceeded"
    ParallelStateAborted = "ParallelStateAborted"
    ParallelStateEntered = "ParallelStateEntered"
    ParallelStateExited = "ParallelStateExited"
    ParallelStateFailed = "ParallelStateFailed"
    ParallelStateStarted = "ParallelStateStarted"
    ParallelStateSucceeded = "ParallelStateSucceeded"
    PassStateEntered = "PassStateEntered"
    PassStateExited = "PassStateExited"
    SucceedStateEntered = "SucceedStateEntered"
    SucceedStateExited = "SucceedStateExited"
    TaskFailed = "TaskFailed"
    TaskScheduled = "TaskScheduled"
    TaskStarted = "TaskStarted"
    TaskStartFailed = "TaskStartFailed"
    TaskStateAborted = "TaskStateAborted"
    TaskStateEntered = "TaskStateEntered"
    TaskStateExited = "TaskStateExited"
    TaskSubmitFailed = "TaskSubmitFailed"
    TaskSubmitted = "TaskSubmitted"
    TaskSucceeded = "TaskSucceeded"
    TaskTimedOut = "TaskTimedOut"
    WaitStateAborted = "WaitStateAborted"
    WaitStateEntered = "WaitStateEntered"
    WaitStateExited = "WaitStateExited"
    MapRunAborted = "MapRunAborted"
    MapRunFailed = "MapRunFailed"
    MapRunStarted = "MapRunStarted"
    MapRunSucceeded = "MapRunSucceeded"
    ExecutionRedriven = "ExecutionRedriven"
    MapRunRedriven = "MapRunRedriven"
    EvaluationFailed = "EvaluationFailed"


class IncludedData(str, Enum):
    ALL_DATA = "ALL_DATA"
    METADATA_ONLY = "METADATA_ONLY"


class InspectionLevel(str, Enum):
    INFO = "INFO"
    DEBUG = "DEBUG"
    TRACE = "TRACE"


class KmsKeyState(str, Enum):
    DISABLED = "DISABLED"
    PENDING_DELETION = "PENDING_DELETION"
    PENDING_IMPORT = "PENDING_IMPORT"
    UNAVAILABLE = "UNAVAILABLE"
    CREATING = "CREATING"


class LogLevel(str, Enum):
    ALL = "ALL"
    ERROR = "ERROR"
    FATAL = "FATAL"
    OFF = "OFF"


class MapRunStatus(str, Enum):
    RUNNING = "RUNNING"
    SUCCEEDED = "SUCCEEDED"
    FAILED = "FAILED"
    ABORTED = "ABORTED"


class StateMachineStatus(str, Enum):
    ACTIVE = "ACTIVE"
    DELETING = "DELETING"


class StateMachineType(str, Enum):
    STANDARD = "STANDARD"
    EXPRESS = "EXPRESS"


class SyncExecutionStatus(str, Enum):
    SUCCEEDED = "SUCCEEDED"
    FAILED = "FAILED"
    TIMED_OUT = "TIMED_OUT"


class TestExecutionStatus(str, Enum):
    SUCCEEDED = "SUCCEEDED"
    FAILED = "FAILED"
    RETRIABLE = "RETRIABLE"
    CAUGHT_ERROR = "CAUGHT_ERROR"


class ValidateStateMachineDefinitionResultCode(str, Enum):
    OK = "OK"
    FAIL = "FAIL"


class ValidateStateMachineDefinitionSeverity(str, Enum):
    ERROR = "ERROR"
    WARNING = "WARNING"


class ValidationExceptionReason(str, Enum):
    API_DOES_NOT_SUPPORT_LABELED_ARNS = "API_DOES_NOT_SUPPORT_LABELED_ARNS"
    MISSING_REQUIRED_PARAMETER = "MISSING_REQUIRED_PARAMETER"
    CANNOT_UPDATE_COMPLETED_MAP_RUN = "CANNOT_UPDATE_COMPLETED_MAP_RUN"
    INVALID_ROUTING_CONFIGURATION = "INVALID_ROUTING_CONFIGURATION"


class ActivityAlreadyExists(ServiceException):
    code: str = "ActivityAlreadyExists"
    sender_fault: bool = False
    status_code: int = 400


class ActivityDoesNotExist(ServiceException):
    code: str = "ActivityDoesNotExist"
    sender_fault: bool = False
    status_code: int = 400


class ActivityLimitExceeded(ServiceException):
    code: str = "ActivityLimitExceeded"
    sender_fault: bool = False
    status_code: int = 400


class ActivityWorkerLimitExceeded(ServiceException):
    code: str = "ActivityWorkerLimitExceeded"
    sender_fault: bool = False
    status_code: int = 400


class ConflictException(ServiceException):
    code: str = "ConflictException"
    sender_fault: bool = False
    status_code: int = 400


class ExecutionAlreadyExists(ServiceException):
    code: str = "ExecutionAlreadyExists"
    sender_fault: bool = False
    status_code: int = 400


class ExecutionDoesNotExist(ServiceException):
    code: str = "ExecutionDoesNotExist"
    sender_fault: bool = False
    status_code: int = 400


class ExecutionLimitExceeded(ServiceException):
    code: str = "ExecutionLimitExceeded"
    sender_fault: bool = False
    status_code: int = 400


class ExecutionNotRedrivable(ServiceException):
    code: str = "ExecutionNotRedrivable"
    sender_fault: bool = False
    status_code: int = 400


class InvalidArn(ServiceException):
    code: str = "InvalidArn"
    sender_fault: bool = False
    status_code: int = 400


class InvalidDefinition(ServiceException):
    code: str = "InvalidDefinition"
    sender_fault: bool = False
    status_code: int = 400


class InvalidEncryptionConfiguration(ServiceException):
    code: str = "InvalidEncryptionConfiguration"
    sender_fault: bool = False
    status_code: int = 400


class InvalidExecutionInput(ServiceException):
    code: str = "InvalidExecutionInput"
    sender_fault: bool = False
    status_code: int = 400


class InvalidLoggingConfiguration(ServiceException):
    code: str = "InvalidLoggingConfiguration"
    sender_fault: bool = False
    status_code: int = 400


class InvalidName(ServiceException):
    code: str = "InvalidName"
    sender_fault: bool = False
    status_code: int = 400


class InvalidOutput(ServiceException):
    code: str = "InvalidOutput"
    sender_fault: bool = False
    status_code: int = 400


class InvalidToken(ServiceException):
    code: str = "InvalidToken"
    exception_type: str = "InvalidToken"
    sender_fault: bool = False
    status_code: int = 400
    message: str = "Invalid Token: 'Invalid token'"

    def __init__(self):
        super().__init__(self.message, self.exception_type, self.status_code)


class InvalidTracingConfiguration(ServiceException):
    code: str = "InvalidTracingConfiguration"
    sender_fault: bool = False
    status_code: int = 400


class KmsAccessDeniedException(ServiceException):
    code: str = "KmsAccessDeniedException"
    sender_fault: bool = False
    status_code: int = 400


class KmsInvalidStateException(ServiceException):
    code: str = "KmsInvalidStateException"
    sender_fault: bool = False
    status_code: int = 400
    kmsKeyState: Optional[KmsKeyState]


class KmsThrottlingException(ServiceException):
    code: str = "KmsThrottlingException"
    sender_fault: bool = False
    status_code: int = 400


class MissingRequiredParameter(ServiceException):
    code: str = "MissingRequiredParameter"
    sender_fault: bool = False
    status_code: int = 400


class ResourceNotFound(ServiceException):
    code: str = "ResourceNotFound"
    sender_fault: bool = False
    status_code: int = 400
    resourceName: Optional[Arn]


class ServiceQuotaExceededException(ServiceException):
    code: str = "ServiceQuotaExceededException"
    sender_fault: bool = False
    status_code: int = 400


class StateMachineAlreadyExists(ServiceException):
    code: str = "StateMachineAlreadyExists"
    sender_fault: bool = False
    status_code: int = 400


class StateMachineDeleting(ServiceException):
    code: str = "StateMachineDeleting"
    sender_fault: bool = False
    status_code: int = 400


class StateMachineDoesNotExist(ServiceException):
    code: str = "StateMachineDoesNotExist"
    sender_fault: bool = False
    status_code: int = 400


class StateMachineLimitExceeded(ServiceException):
    code: str = "StateMachineLimitExceeded"
    sender_fault: bool = False
    status_code: int = 400


class StateMachineTypeNotSupported(ServiceException):
    code: str = "StateMachineTypeNotSupported"
    sender_fault: bool = False
    status_code: int = 400


class TaskDoesNotExist(ServiceException):
    code: str = "TaskDoesNotExist"
    sender_fault: bool = False
    status_code: int = 400


class TaskTimedOut(ServiceException):
    code: str = "TaskTimedOut"
    sender_fault: bool = False
    status_code: int = 400


class TooManyTags(ServiceException):
    code: str = "TooManyTags"
    sender_fault: bool = False
    status_code: int = 400
    resourceName: Optional[Arn]


class ValidationException(ServiceException):
    code: str = "ValidationException"
    sender_fault: bool = False
    status_code: int = 400
    reason: Optional[ValidationExceptionReason]


class ActivityFailedEventDetails(TypedDict, total=False):
    error: Optional[SensitiveError]
    cause: Optional[SensitiveCause]


Timestamp = datetime


class ActivityListItem(TypedDict, total=False):
    activityArn: Arn
    name: Name
    creationDate: Timestamp


ActivityList = List[ActivityListItem]


class ActivityScheduleFailedEventDetails(TypedDict, total=False):
    error: Optional[SensitiveError]
    cause: Optional[SensitiveCause]


TimeoutInSeconds = int


class HistoryEventExecutionDataDetails(TypedDict, total=False):
    truncated: Optional[truncated]


class ActivityScheduledEventDetails(TypedDict, total=False):
    resource: Arn
    input: Optional[SensitiveData]
    inputDetails: Optional[HistoryEventExecutionDataDetails]
    timeoutInSeconds: Optional[TimeoutInSeconds]
    heartbeatInSeconds: Optional[TimeoutInSeconds]


class ActivityStartedEventDetails(TypedDict, total=False):
    workerName: Optional[Identity]


class ActivitySucceededEventDetails(TypedDict, total=False):
    output: Optional[SensitiveData]
    outputDetails: Optional[HistoryEventExecutionDataDetails]


class ActivityTimedOutEventDetails(TypedDict, total=False):
    error: Optional[SensitiveError]
    cause: Optional[SensitiveCause]


AssignedVariables = Dict[VariableName, VariableValue]


class AssignedVariablesDetails(TypedDict, total=False):
    truncated: Optional[truncated]


BilledDuration = int
BilledMemoryUsed = int


class BillingDetails(TypedDict, total=False):
    billedMemoryUsedInMB: Optional[BilledMemoryUsed]
    billedDurationInMilliseconds: Optional[BilledDuration]


class CloudWatchEventsExecutionDataDetails(TypedDict, total=False):
    included: Optional[includedDetails]


class CloudWatchLogsLogGroup(TypedDict, total=False):
    logGroupArn: Optional[Arn]


EncryptionConfiguration = TypedDict(
    "EncryptionConfiguration",
    {
        "kmsKeyId": Optional[KmsKeyId],
        "kmsDataKeyReusePeriodSeconds": Optional[KmsDataKeyReusePeriodSeconds],
        "type": EncryptionType,
    },
    total=False,
)


class Tag(TypedDict, total=False):
    key: Optional[TagKey]
    value: Optional[TagValue]


TagList = List[Tag]


class CreateActivityOutput(TypedDict, total=False):
    activityArn: Arn
    creationDate: Timestamp


class RoutingConfigurationListItem(TypedDict, total=False):
    stateMachineVersionArn: Arn
    weight: VersionWeight


RoutingConfigurationList = List[RoutingConfigurationListItem]


class CreateStateMachineAliasOutput(TypedDict, total=False):
    stateMachineAliasArn: Arn
    creationDate: Timestamp


class TracingConfiguration(TypedDict, total=False):
    enabled: Optional[Enabled]


class LogDestination(TypedDict, total=False):
    cloudWatchLogsLogGroup: Optional[CloudWatchLogsLogGroup]


LogDestinationList = List[LogDestination]


class LoggingConfiguration(TypedDict, total=False):
    level: Optional[LogLevel]
    includeExecutionData: Optional[IncludeExecutionData]
    destinations: Optional[LogDestinationList]


CreateStateMachineInput = TypedDict(
    "CreateStateMachineInput",
    {
        "name": Name,
        "definition": Definition,
        "roleArn": Arn,
        "type": Optional[StateMachineType],
        "loggingConfiguration": Optional[LoggingConfiguration],
        "tags": Optional[TagList],
        "tracingConfiguration": Optional[TracingConfiguration],
        "publish": Optional[Publish],
        "versionDescription": Optional[VersionDescription],
        "encryptionConfiguration": Optional[EncryptionConfiguration],
    },
    total=False,
)


class CreateStateMachineOutput(TypedDict, total=False):
    stateMachineArn: Arn
    creationDate: Timestamp
    stateMachineVersionArn: Optional[Arn]


class DeleteActivityOutput(TypedDict, total=False):
    pass


class DeleteStateMachineAliasOutput(TypedDict, total=False):
    pass


class DeleteStateMachineOutput(TypedDict, total=False):
    pass


class DeleteStateMachineVersionOutput(TypedDict, total=False):
    pass


class DescribeActivityOutput(TypedDict, total=False):
    activityArn: Arn
    name: Name
    creationDate: Timestamp
    encryptionConfiguration: Optional[EncryptionConfiguration]


class DescribeExecutionOutput(TypedDict, total=False):
    executionArn: Arn
    stateMachineArn: Arn
    name: Optional[Name]
    status: ExecutionStatus
    startDate: Timestamp
    stopDate: Optional[Timestamp]
    input: Optional[SensitiveData]
    inputDetails: Optional[CloudWatchEventsExecutionDataDetails]
    output: Optional[SensitiveData]
    outputDetails: Optional[CloudWatchEventsExecutionDataDetails]
    traceHeader: Optional[TraceHeader]
    mapRunArn: Optional[LongArn]
    error: Optional[SensitiveError]
    cause: Optional[SensitiveCause]
    stateMachineVersionArn: Optional[Arn]
    stateMachineAliasArn: Optional[Arn]
    redriveCount: Optional[RedriveCount]
    redriveDate: Optional[Timestamp]
    redriveStatus: Optional[ExecutionRedriveStatus]
    redriveStatusReason: Optional[SensitiveData]


LongObject = int
UnsignedLong = int


class MapRunExecutionCounts(TypedDict, total=False):
    pending: UnsignedLong
    running: UnsignedLong
    succeeded: UnsignedLong
    failed: UnsignedLong
    timedOut: UnsignedLong
    aborted: UnsignedLong
    total: UnsignedLong
    resultsWritten: UnsignedLong
    failuresNotRedrivable: Optional[LongObject]
    pendingRedrive: Optional[LongObject]


class MapRunItemCounts(TypedDict, total=False):
    pending: UnsignedLong
    running: UnsignedLong
    succeeded: UnsignedLong
    failed: UnsignedLong
    timedOut: UnsignedLong
    aborted: UnsignedLong
    total: UnsignedLong
    resultsWritten: UnsignedLong
    failuresNotRedrivable: Optional[LongObject]
    pendingRedrive: Optional[LongObject]


ToleratedFailureCount = int


class DescribeMapRunOutput(TypedDict, total=False):
    mapRunArn: LongArn
    executionArn: Arn
    status: MapRunStatus
    startDate: Timestamp
    stopDate: Optional[Timestamp]
    maxConcurrency: MaxConcurrency
    toleratedFailurePercentage: ToleratedFailurePercentage
    toleratedFailureCount: ToleratedFailureCount
    itemCounts: MapRunItemCounts
    executionCounts: MapRunExecutionCounts
    redriveCount: Optional[RedriveCount]
    redriveDate: Optional[Timestamp]


class DescribeStateMachineAliasOutput(TypedDict, total=False):
    stateMachineAliasArn: Optional[Arn]
    name: Optional[Name]
    description: Optional[AliasDescription]
    routingConfiguration: Optional[RoutingConfigurationList]
    creationDate: Optional[Timestamp]
    updateDate: Optional[Timestamp]


VariableNameList = List[VariableName]
VariableReferences = Dict[StateName, VariableNameList]


class DescribeStateMachineForExecutionOutput(TypedDict, total=False):
    stateMachineArn: Arn
    name: Name
    definition: Definition
    roleArn: Arn
    updateDate: Timestamp
    loggingConfiguration: Optional[LoggingConfiguration]
    tracingConfiguration: Optional[TracingConfiguration]
    mapRunArn: Optional[LongArn]
    label: Optional[MapRunLabel]
    revisionId: Optional[RevisionId]
    encryptionConfiguration: Optional[EncryptionConfiguration]
    variableReferences: Optional[VariableReferences]


DescribeStateMachineOutput = TypedDict(
    "DescribeStateMachineOutput",
    {
        "stateMachineArn": Arn,
        "name": Name,
        "status": Optional[StateMachineStatus],
        "definition": Definition,
        "roleArn": Arn,
        "type": StateMachineType,
        "creationDate": Timestamp,
        "loggingConfiguration": Optional[LoggingConfiguration],
        "tracingConfiguration": Optional[TracingConfiguration],
        "label": Optional[MapRunLabel],
        "revisionId": Optional[RevisionId],
        "description": Optional[VersionDescription],
        "encryptionConfiguration": Optional[EncryptionConfiguration],
        "variableReferences": Optional[VariableReferences],
    },
    total=False,
)


class EvaluationFailedEventDetails(TypedDict, total=False):
    error: Optional[SensitiveError]
    cause: Optional[SensitiveCause]
    location: Optional[EvaluationFailureLocation]
    state: StateName


EventId = int


class ExecutionAbortedEventDetails(TypedDict, total=False):
    error: Optional[SensitiveError]
    cause: Optional[SensitiveCause]


class ExecutionFailedEventDetails(TypedDict, total=False):
    error: Optional[SensitiveError]
    cause: Optional[SensitiveCause]


class ExecutionListItem(TypedDict, total=False):
    executionArn: Arn
    stateMachineArn: Arn
    name: Name
    status: ExecutionStatus
    startDate: Timestamp
    stopDate: Optional[Timestamp]
    mapRunArn: Optional[LongArn]
    itemCount: Optional[UnsignedInteger]
    stateMachineVersionArn: Optional[Arn]
    stateMachineAliasArn: Optional[Arn]
    redriveCount: Optional[RedriveCount]
    redriveDate: Optional[Timestamp]


ExecutionList = List[ExecutionListItem]


class ExecutionRedrivenEventDetails(TypedDict, total=False):
    redriveCount: Optional[RedriveCount]


class ExecutionStartedEventDetails(TypedDict, total=False):
    input: Optional[SensitiveData]
    inputDetails: Optional[HistoryEventExecutionDataDetails]
    roleArn: Optional[Arn]
    stateMachineAliasArn: Optional[Arn]
    stateMachineVersionArn: Optional[Arn]


class ExecutionSucceededEventDetails(TypedDict, total=False):
    output: Optional[SensitiveData]
    outputDetails: Optional[HistoryEventExecutionDataDetails]


class ExecutionTimedOutEventDetails(TypedDict, total=False):
    error: Optional[SensitiveError]
    cause: Optional[SensitiveCause]


class GetActivityTaskOutput(TypedDict, total=False):
    taskToken: Optional[TaskToken]
    input: Optional[SensitiveDataJobInput]


class MapRunRedrivenEventDetails(TypedDict, total=False):
    mapRunArn: Optional[LongArn]
    redriveCount: Optional[RedriveCount]


class MapRunFailedEventDetails(TypedDict, total=False):
    error: Optional[SensitiveError]
    cause: Optional[SensitiveCause]


class MapRunStartedEventDetails(TypedDict, total=False):
    mapRunArn: Optional[LongArn]


class StateExitedEventDetails(TypedDict, total=False):
    name: Name
    output: Optional[SensitiveData]
    outputDetails: Optional[HistoryEventExecutionDataDetails]
    assignedVariables: Optional[AssignedVariables]
    assignedVariablesDetails: Optional[AssignedVariablesDetails]


class StateEnteredEventDetails(TypedDict, total=False):
    name: Name
    input: Optional[SensitiveData]
    inputDetails: Optional[HistoryEventExecutionDataDetails]


class LambdaFunctionTimedOutEventDetails(TypedDict, total=False):
    error: Optional[SensitiveError]
    cause: Optional[SensitiveCause]


class LambdaFunctionSucceededEventDetails(TypedDict, total=False):
    output: Optional[SensitiveData]
    outputDetails: Optional[HistoryEventExecutionDataDetails]


class LambdaFunctionStartFailedEventDetails(TypedDict, total=False):
    error: Optional[SensitiveError]
    cause: Optional[SensitiveCause]


class TaskCredentials(TypedDict, total=False):
    roleArn: Optional[LongArn]


class LambdaFunctionScheduledEventDetails(TypedDict, total=False):
    resource: Arn
    input: Optional[SensitiveData]
    inputDetails: Optional[HistoryEventExecutionDataDetails]
    timeoutInSeconds: Optional[TimeoutInSeconds]
    taskCredentials: Optional[TaskCredentials]


class LambdaFunctionScheduleFailedEventDetails(TypedDict, total=False):
    error: Optional[SensitiveError]
    cause: Optional[SensitiveCause]


class LambdaFunctionFailedEventDetails(TypedDict, total=False):
    error: Optional[SensitiveError]
    cause: Optional[SensitiveCause]


class MapIterationEventDetails(TypedDict, total=False):
    name: Optional[Name]
    index: Optional[UnsignedInteger]


class MapStateStartedEventDetails(TypedDict, total=False):
    length: Optional[UnsignedInteger]


class TaskTimedOutEventDetails(TypedDict, total=False):
    resourceType: Name
    resource: Name
    error: Optional[SensitiveError]
    cause: Optional[SensitiveCause]


class TaskSucceededEventDetails(TypedDict, total=False):
    resourceType: Name
    resource: Name
    output: Optional[SensitiveData]
    outputDetails: Optional[HistoryEventExecutionDataDetails]


class TaskSubmittedEventDetails(TypedDict, total=False):
    resourceType: Name
    resource: Name
    output: Optional[SensitiveData]
    outputDetails: Optional[HistoryEventExecutionDataDetails]


class TaskSubmitFailedEventDetails(TypedDict, total=False):
    resourceType: Name
    resource: Name
    error: Optional[SensitiveError]
    cause: Optional[SensitiveCause]


class TaskStartedEventDetails(TypedDict, total=False):
    resourceType: Name
    resource: Name


class TaskStartFailedEventDetails(TypedDict, total=False):
    resourceType: Name
    resource: Name
    error: Optional[SensitiveError]
    cause: Optional[SensitiveCause]


class TaskScheduledEventDetails(TypedDict, total=False):
    resourceType: Name
    resource: Name
    region: Name
    parameters: ConnectorParameters
    timeoutInSeconds: Optional[TimeoutInSeconds]
    heartbeatInSeconds: Optional[TimeoutInSeconds]
    taskCredentials: Optional[TaskCredentials]


class TaskFailedEventDetails(TypedDict, total=False):
    resourceType: Name
    resource: Name
    error: Optional[SensitiveError]
    cause: Optional[SensitiveCause]


HistoryEvent = TypedDict(
    "HistoryEvent",
    {
        "timestamp": Timestamp,
        "type": HistoryEventType,
        "id": EventId,
        "previousEventId": Optional[EventId],
        "activityFailedEventDetails": Optional[ActivityFailedEventDetails],
        "activityScheduleFailedEventDetails": Optional[
            ActivityScheduleFailedEventDetails
        ],
        "activityScheduledEventDetails": Optional[ActivityScheduledEventDetails],
        "activityStartedEventDetails": Optional[ActivityStartedEventDetails],
        "activitySucceededEventDetails": Optional[ActivitySucceededEventDetails],
        "activityTimedOutEventDetails": Optional[ActivityTimedOutEventDetails],
        "taskFailedEventDetails": Optional[TaskFailedEventDetails],
        "taskScheduledEventDetails": Optional[TaskScheduledEventDetails],
        "taskStartFailedEventDetails": Optional[TaskStartFailedEventDetails],
        "taskStartedEventDetails": Optional[TaskStartedEventDetails],
        "taskSubmitFailedEventDetails": Optional[TaskSubmitFailedEventDetails],
        "taskSubmittedEventDetails": Optional[TaskSubmittedEventDetails],
        "taskSucceededEventDetails": Optional[TaskSucceededEventDetails],
        "taskTimedOutEventDetails": Optional[TaskTimedOutEventDetails],
        "executionFailedEventDetails": Optional[ExecutionFailedEventDetails],
        "executionStartedEventDetails": Optional[ExecutionStartedEventDetails],
        "executionSucceededEventDetails": Optional[ExecutionSucceededEventDetails],
        "executionAbortedEventDetails": Optional[ExecutionAbortedEventDetails],
        "executionTimedOutEventDetails": Optional[ExecutionTimedOutEventDetails],
        "executionRedrivenEventDetails": Optional[ExecutionRedrivenEventDetails],
        "mapStateStartedEventDetails": Optional[MapStateStartedEventDetails],
        "mapIterationStartedEventDetails": Optional[MapIterationEventDetails],
        "mapIterationSucceededEventDetails": Optional[MapIterationEventDetails],
        "mapIterationFailedEventDetails": Optional[MapIterationEventDetails],
        "mapIterationAbortedEventDetails": Optional[MapIterationEventDetails],
        "lambdaFunctionFailedEventDetails": Optional[LambdaFunctionFailedEventDetails],
        "lambdaFunctionScheduleFailedEventDetails": Optional[
            LambdaFunctionScheduleFailedEventDetails
        ],
        "lambdaFunctionScheduledEventDetails": Optional[
            LambdaFunctionScheduledEventDetails
        ],
        "lambdaFunctionStartFailedEventDetails": Optional[
            LambdaFunctionStartFailedEventDetails
        ],
        "lambdaFunctionSucceededEventDetails": Optional[
            LambdaFunctionSucceededEventDetails
        ],
        "lambdaFunctionTimedOutEventDetails": Optional[
            LambdaFunctionTimedOutEventDetails
        ],
        "stateEnteredEventDetails": Optional[StateEnteredEventDetails],
        "stateExitedEventDetails": Optional[StateExitedEventDetails],
        "mapRunStartedEventDetails": Optional[MapRunStartedEventDetails],
        "mapRunFailedEventDetails": Optional[MapRunFailedEventDetails],
        "mapRunRedrivenEventDetails": Optional[MapRunRedrivenEventDetails],
        "evaluationFailedEventDetails": Optional[EvaluationFailedEventDetails],
    },
    total=False,
)
HistoryEventList = List[HistoryEvent]


class GetExecutionHistoryOutput(TypedDict, total=False):
    events: HistoryEventList
    nextToken: Optional[PageToken]


class InspectionDataResponse(TypedDict, total=False):
    protocol: Optional[HTTPProtocol]
    statusCode: Optional[HTTPStatusCode]
    statusMessage: Optional[HTTPStatusMessage]
    headers: Optional[HTTPHeaders]
    body: Optional[HTTPBody]


class InspectionDataRequest(TypedDict, total=False):
    protocol: Optional[HTTPProtocol]
    method: Optional[HTTPMethod]
    url: Optional[URL]
    headers: Optional[HTTPHeaders]
    body: Optional[HTTPBody]


class InspectionData(TypedDict, total=False):
    input: Optional[SensitiveData]
    afterArguments: Optional[SensitiveData]
    afterInputPath: Optional[SensitiveData]
    afterParameters: Optional[SensitiveData]
    result: Optional[SensitiveData]
    afterResultSelector: Optional[SensitiveData]
    afterResultPath: Optional[SensitiveData]
    request: Optional[InspectionDataRequest]
    response: Optional[InspectionDataResponse]
    variables: Optional[SensitiveData]


class MapRunListItem(TypedDict, total=False):
    executionArn: Arn
    mapRunArn: LongArn
    stateMachineArn: Arn
    startDate: Timestamp
    stopDate: Optional[Timestamp]


MapRunList = List[MapRunListItem]


class ListMapRunsOutput(TypedDict, total=False):
    mapRuns: MapRunList
    nextToken: Optional[PageToken]


class StateMachineAliasListItem(TypedDict, total=False):
    stateMachineAliasArn: LongArn
    creationDate: Timestamp


StateMachineAliasList = List[StateMachineAliasListItem]


class ListStateMachineAliasesOutput(TypedDict, total=False):
    stateMachineAliases: StateMachineAliasList
    nextToken: Optional[PageToken]


class StateMachineVersionListItem(TypedDict, total=False):
    stateMachineVersionArn: LongArn
    creationDate: Timestamp


StateMachineVersionList = List[StateMachineVersionListItem]


class ListStateMachineVersionsOutput(TypedDict, total=False):
    stateMachineVersions: StateMachineVersionList
    nextToken: Optional[PageToken]


StateMachineListItem = TypedDict(
    "StateMachineListItem",
    {
        "stateMachineArn": Arn,
        "name": Name,
        "type": StateMachineType,
        "creationDate": Timestamp,
    },
    total=False,
)
StateMachineList = List[StateMachineListItem]


class ListStateMachinesOutput(TypedDict, total=False):
    stateMachines: StateMachineList
    nextToken: Optional[PageToken]


class ListTagsForResourceOutput(TypedDict, total=False):
    tags: Optional[TagList]


class PublishStateMachineVersionOutput(TypedDict, total=False):
    creationDate: Timestamp
    stateMachineVersionArn: Arn


class RedriveExecutionOutput(TypedDict, total=False):
    redriveDate: Timestamp


class SendTaskFailureOutput(TypedDict, total=False):
    pass


class SendTaskHeartbeatOutput(TypedDict, total=False):
    pass


class SendTaskSuccessOutput(TypedDict, total=False):
    pass


class StartExecutionOutput(TypedDict, total=False):
    executionArn: Arn
    startDate: Timestamp


class StartSyncExecutionOutput(TypedDict, total=False):
    executionArn: Arn
    stateMachineArn: Optional[Arn]
    name: Optional[Name]
    startDate: Timestamp
    stopDate: Timestamp
    status: SyncExecutionStatus
    error: Optional[SensitiveError]
    cause: Optional[SensitiveCause]
    input: Optional[SensitiveData]
    inputDetails: Optional[CloudWatchEventsExecutionDataDetails]
    output: Optional[SensitiveData]
    outputDetails: Optional[CloudWatchEventsExecutionDataDetails]
    traceHeader: Optional[TraceHeader]
    billingDetails: Optional[BillingDetails]


class StopExecutionOutput(TypedDict, total=False):
    stopDate: Timestamp


TagKeyList = List[TagKey]


class TagResourceOutput(TypedDict, total=False):
    pass


class TestStateOutput(TypedDict, total=False):
    output: Optional[SensitiveData]
    error: Optional[SensitiveError]
    cause: Optional[SensitiveCause]
    inspectionData: Optional[InspectionData]
    nextState: Optional[StateName]
    status: Optional[TestExecutionStatus]


class UntagResourceOutput(TypedDict, total=False):
    pass


class UpdateMapRunOutput(TypedDict, total=False):
    pass


class UpdateStateMachineAliasOutput(TypedDict, total=False):
    updateDate: Timestamp


class UpdateStateMachineOutput(TypedDict, total=False):
    updateDate: Timestamp
    revisionId: Optional[RevisionId]
    stateMachineVersionArn: Optional[Arn]


class ValidateStateMachineDefinitionDiagnostic(TypedDict, total=False):
    severity: ValidateStateMachineDefinitionSeverity
    code: ValidateStateMachineDefinitionCode
    message: ValidateStateMachineDefinitionMessage
    location: Optional[ValidateStateMachineDefinitionLocation]


ValidateStateMachineDefinitionDiagnosticList = List[
    ValidateStateMachineDefinitionDiagnostic
]
ValidateStateMachineDefinitionInput = TypedDict(
    "ValidateStateMachineDefinitionInput",
    {
        "definition": Definition,
        "type": Optional[StateMachineType],
        "severity": Optional[ValidateStateMachineDefinitionSeverity],
        "maxResults": Optional[ValidateStateMachineDefinitionMaxResult],
    },
    total=False,
)


class ValidateStateMachineDefinitionOutput(TypedDict, total=False):
    result: ValidateStateMachineDefinitionResultCode
    diagnostics: ValidateStateMachineDefinitionDiagnosticList
    truncated: Optional[ValidateStateMachineDefinitionTruncated]


class InvocationResponse(TypedDict, total=False):
    Payload: Optional[Union[bytes, IO[bytes], Iterable[bytes]]]
    StatusCode: Optional[int]
    FunctionError: Optional[str]
    LogResult: Optional[str]
    ExecutedVersion: Optional[str]
