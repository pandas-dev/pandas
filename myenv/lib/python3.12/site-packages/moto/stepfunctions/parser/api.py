from datetime import datetime
from typing import List, Optional, TypedDict

from ..exceptions import AWSError as ServiceException

AliasDescription = str
CharacterRestrictedName = str
ClientToken = str
ConnectorParameters = str
Definition = str
Enabled = bool
ErrorMessage = str
HTTPBody = str
HTTPHeaders = str
HTTPMethod = str
HTTPProtocol = str
HTTPStatusCode = str
HTTPStatusMessage = str
Identity = str
IncludeExecutionData = bool
IncludeExecutionDataGetExecutionHistory = bool
ListExecutionsPageToken = str
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
VersionDescription = str
VersionWeight = int
includedDetails = bool
truncated = bool


class ExecutionRedriveFilter(str):
    REDRIVEN = "REDRIVEN"
    NOT_REDRIVEN = "NOT_REDRIVEN"


class ExecutionRedriveStatus(str):
    REDRIVABLE = "REDRIVABLE"
    NOT_REDRIVABLE = "NOT_REDRIVABLE"
    REDRIVABLE_BY_MAP_RUN = "REDRIVABLE_BY_MAP_RUN"


class ExecutionStatus(str):
    RUNNING = "RUNNING"
    SUCCEEDED = "SUCCEEDED"
    FAILED = "FAILED"
    TIMED_OUT = "TIMED_OUT"
    ABORTED = "ABORTED"
    PENDING_REDRIVE = "PENDING_REDRIVE"


class HistoryEventType(str):
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


class InspectionLevel(str):
    INFO = "INFO"
    DEBUG = "DEBUG"
    TRACE = "TRACE"


class LogLevel(str):
    ALL = "ALL"
    ERROR = "ERROR"
    FATAL = "FATAL"
    OFF = "OFF"


class MapRunStatus(str):
    RUNNING = "RUNNING"
    SUCCEEDED = "SUCCEEDED"
    FAILED = "FAILED"
    ABORTED = "ABORTED"


class StateMachineStatus(str):
    ACTIVE = "ACTIVE"
    DELETING = "DELETING"


class StateMachineType(str):
    STANDARD = "STANDARD"
    EXPRESS = "EXPRESS"


class SyncExecutionStatus(str):
    SUCCEEDED = "SUCCEEDED"
    FAILED = "FAILED"
    TIMED_OUT = "TIMED_OUT"


class TestExecutionStatus(str):
    SUCCEEDED = "SUCCEEDED"
    FAILED = "FAILED"
    RETRIABLE = "RETRIABLE"
    CAUGHT_ERROR = "CAUGHT_ERROR"


class ValidationExceptionReason(str):
    API_DOES_NOT_SUPPORT_LABELED_ARNS = "API_DOES_NOT_SUPPORT_LABELED_ARNS"
    MISSING_REQUIRED_PARAMETER = "MISSING_REQUIRED_PARAMETER"
    CANNOT_UPDATE_COMPLETED_MAP_RUN = "CANNOT_UPDATE_COMPLETED_MAP_RUN"
    INVALID_ROUTING_CONFIGURATION = "INVALID_ROUTING_CONFIGURATION"


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
    exception_type: str = "UnrecognizedClientException"
    sender_fault: bool = False
    status_code: int = 400
    message: str = "The security token included in the request is invalid."

    def __init__(self):
        super().__init__(self.message, self.exception_type, self.status_code)


class InvalidTracingConfiguration(ServiceException):
    code: str = "InvalidTracingConfiguration"
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
    resourceName: Optional[str]


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
    activityArn: str
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
    resource: str
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


BilledDuration = int
BilledMemoryUsed = int


class BillingDetails(TypedDict, total=False):
    billedMemoryUsedInMB: Optional[BilledMemoryUsed]
    billedDurationInMilliseconds: Optional[BilledDuration]


class CloudWatchEventsExecutionDataDetails(TypedDict, total=False):
    included: Optional[includedDetails]


class CloudWatchLogsLogGroup(TypedDict, total=False):
    logGroupArn: Optional[str]


class Tag(TypedDict, total=False):
    key: Optional[TagKey]
    value: Optional[TagValue]


TagList = List[Tag]


class CreateActivityOutput(TypedDict, total=False):
    activityArn: str
    creationDate: Timestamp


class RoutingConfigurationListItem(TypedDict, total=False):
    stateMachineVersionArn: str
    weight: VersionWeight


RoutingConfigurationList = List[RoutingConfigurationListItem]


class CreateStateMachineAliasOutput(TypedDict, total=False):
    stateMachineAliasArn: str
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
        "roleArn": str,
        "type": Optional[StateMachineType],
        "loggingConfiguration": Optional[LoggingConfiguration],
        "tags": Optional[TagList],
        "tracingConfiguration": Optional[TracingConfiguration],
        "publish": Optional[Publish],
        "versionDescription": Optional[VersionDescription],
    },
    total=False,
)


class CreateStateMachineOutput(TypedDict, total=False):
    stateMachineArn: str
    creationDate: Timestamp
    stateMachineVersionArn: Optional[str]


class DeleteActivityOutput(TypedDict, total=False):
    pass


class DeleteStateMachineAliasOutput(TypedDict, total=False):
    pass


class DeleteStateMachineOutput(TypedDict, total=False):
    pass


class DeleteStateMachineVersionOutput(TypedDict, total=False):
    pass


class DescribeActivityOutput(TypedDict, total=False):
    activityArn: str
    name: Name
    creationDate: Timestamp


class DescribeExecutionOutput(TypedDict, total=False):
    executionArn: str
    stateMachineArn: str
    name: Optional[Name]
    status: ExecutionStatus
    startDate: Timestamp
    stopDate: Optional[Timestamp]
    input: Optional[SensitiveData]
    inputDetails: Optional[CloudWatchEventsExecutionDataDetails]
    output: Optional[SensitiveData]
    outputDetails: Optional[CloudWatchEventsExecutionDataDetails]
    traceHeader: Optional[TraceHeader]
    mapRunArn: Optional[str]
    error: Optional[SensitiveError]
    cause: Optional[SensitiveCause]
    stateMachineVersionArn: Optional[str]
    stateMachineAliasArn: Optional[str]
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
    mapRunArn: str
    executionArn: str
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
    stateMachineAliasArn: Optional[str]
    name: Optional[Name]
    description: Optional[AliasDescription]
    routingConfiguration: Optional[RoutingConfigurationList]
    creationDate: Optional[Timestamp]
    updateDate: Optional[Timestamp]


class DescribeStateMachineForExecutionOutput(TypedDict, total=False):
    stateMachineArn: str
    name: Name
    definition: Definition
    roleArn: str
    updateDate: Timestamp
    loggingConfiguration: Optional[LoggingConfiguration]
    tracingConfiguration: Optional[TracingConfiguration]
    mapRunArn: Optional[str]
    label: Optional[MapRunLabel]
    revisionId: Optional[RevisionId]


DescribeStateMachineOutput = TypedDict(
    "DescribeStateMachineOutput",
    {
        "stateMachineArn": str,
        "name": Name,
        "status": Optional[StateMachineStatus],
        "definition": Definition,
        "roleArn": str,
        "type": StateMachineType,
        "creationDate": Timestamp,
        "loggingConfiguration": Optional[LoggingConfiguration],
        "tracingConfiguration": Optional[TracingConfiguration],
        "label": Optional[MapRunLabel],
        "revisionId": Optional[RevisionId],
        "description": Optional[VersionDescription],
    },
    total=False,
)
EventId = int


class ExecutionAbortedEventDetails(TypedDict, total=False):
    error: Optional[SensitiveError]
    cause: Optional[SensitiveCause]


class ExecutionFailedEventDetails(TypedDict, total=False):
    error: Optional[SensitiveError]
    cause: Optional[SensitiveCause]


class ExecutionListItem(TypedDict, total=False):
    executionArn: str
    stateMachineArn: str
    name: Name
    status: ExecutionStatus
    startDate: Timestamp
    stopDate: Optional[Timestamp]
    mapRunArn: Optional[str]
    itemCount: Optional[UnsignedInteger]
    stateMachineVersionArn: Optional[str]
    stateMachineAliasArn: Optional[str]
    redriveCount: Optional[RedriveCount]
    redriveDate: Optional[Timestamp]


ExecutionList = List[ExecutionListItem]


class ExecutionRedrivenEventDetails(TypedDict, total=False):
    redriveCount: Optional[RedriveCount]


class ExecutionStartedEventDetails(TypedDict, total=False):
    input: Optional[SensitiveData]
    inputDetails: Optional[HistoryEventExecutionDataDetails]
    roleArn: Optional[str]
    stateMachineAliasArn: Optional[str]
    stateMachineVersionArn: Optional[str]


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
    mapRunArn: Optional[str]
    redriveCount: Optional[RedriveCount]


class MapRunFailedEventDetails(TypedDict, total=False):
    error: Optional[SensitiveError]
    cause: Optional[SensitiveCause]


class MapRunStartedEventDetails(TypedDict, total=False):
    mapRunArn: Optional[str]


class StateExitedEventDetails(TypedDict, total=False):
    name: Name
    output: Optional[SensitiveData]
    outputDetails: Optional[HistoryEventExecutionDataDetails]


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
    roleArn: Optional[str]


class LambdaFunctionScheduledEventDetails(TypedDict, total=False):
    resource: str
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
    afterInputPath: Optional[SensitiveData]
    afterParameters: Optional[SensitiveData]
    result: Optional[SensitiveData]
    afterResultSelector: Optional[SensitiveData]
    afterResultPath: Optional[SensitiveData]
    request: Optional[InspectionDataRequest]
    response: Optional[InspectionDataResponse]


class ListActivitiesOutput(TypedDict, total=False):
    activities: ActivityList
    nextToken: Optional[PageToken]


class ListExecutionsOutput(TypedDict, total=False):
    executions: ExecutionList
    nextToken: Optional[ListExecutionsPageToken]


class MapRunListItem(TypedDict, total=False):
    executionArn: str
    mapRunArn: str
    stateMachineArn: str
    startDate: Timestamp
    stopDate: Optional[Timestamp]


MapRunList = List[MapRunListItem]


class ListMapRunsOutput(TypedDict, total=False):
    mapRuns: MapRunList
    nextToken: Optional[PageToken]


class StateMachineAliasListItem(TypedDict, total=False):
    stateMachineAliasArn: str
    creationDate: Timestamp


StateMachineAliasList = List[StateMachineAliasListItem]


class ListStateMachineAliasesOutput(TypedDict, total=False):
    stateMachineAliases: StateMachineAliasList
    nextToken: Optional[PageToken]


class StateMachineVersionListItem(TypedDict, total=False):
    stateMachineVersionArn: str
    creationDate: Timestamp


StateMachineVersionList = List[StateMachineVersionListItem]


class ListStateMachineVersionsOutput(TypedDict, total=False):
    stateMachineVersions: StateMachineVersionList
    nextToken: Optional[PageToken]


StateMachineListItem = TypedDict(
    "StateMachineListItem",
    {
        "stateMachineArn": str,
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
    stateMachineVersionArn: str


class RedriveExecutionOutput(TypedDict, total=False):
    redriveDate: Timestamp


class SendTaskFailureOutput(TypedDict, total=False):
    pass


class SendTaskHeartbeatOutput(TypedDict, total=False):
    pass


class SendTaskSuccessOutput(TypedDict, total=False):
    pass


class StartExecutionOutput(TypedDict, total=False):
    executionArn: str
    startDate: Timestamp


class StartSyncExecutionOutput(TypedDict, total=False):
    executionArn: str
    stateMachineArn: Optional[str]
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
    stateMachineVersionArn: Optional[str]
