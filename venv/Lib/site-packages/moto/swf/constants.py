# List decision fields and if they're required or not
#
# See http://docs.aws.amazon.com/amazonswf/latest/apireference/API_RespondDecisionTaskCompleted.html
# and subsequent docs for each decision type.
DECISIONS_FIELDS = {
    "cancelTimerDecisionAttributes": {"timerId": {"type": "string", "required": True}},
    "cancelWorkflowExecutionDecisionAttributes": {
        "details": {"type": "string", "required": False}
    },
    "completeWorkflowExecutionDecisionAttributes": {
        "result": {"type": "string", "required": False}
    },
    "continueAsNewWorkflowExecutionDecisionAttributes": {
        "childPolicy": {"type": "string", "required": False},
        "executionStartToCloseTimeout": {"type": "string", "required": False},
        "input": {"type": "string", "required": False},
        "lambdaRole": {"type": "string", "required": False},
        "tagList": {"type": "string", "array": True, "required": False},
        "taskList": {"type": "TaskList", "required": False},
        "taskPriority": {"type": "string", "required": False},
        "taskStartToCloseTimeout": {"type": "string", "required": False},
        "workflowTypeVersion": {"type": "string", "required": False},
    },
    "failWorkflowExecutionDecisionAttributes": {
        "details": {"type": "string", "required": False},
        "reason": {"type": "string", "required": False},
    },
    "recordMarkerDecisionAttributes": {
        "details": {"type": "string", "required": False},
        "markerName": {"type": "string", "required": True},
    },
    "requestCancelActivityTaskDecisionAttributes": {
        "activityId": {"type": "string", "required": True}
    },
    "requestCancelExternalWorkflowExecutionDecisionAttributes": {
        "control": {"type": "string", "required": False},
        "runId": {"type": "string", "required": False},
        "workflowId": {"type": "string", "required": True},
    },
    "scheduleActivityTaskDecisionAttributes": {
        "activityId": {"type": "string", "required": True},
        "activityType": {"type": "ActivityType", "required": True},
        "control": {"type": "string", "required": False},
        "heartbeatTimeout": {"type": "string", "required": False},
        "input": {"type": "string", "required": False},
        "scheduleToCloseTimeout": {"type": "string", "required": False},
        "scheduleToStartTimeout": {"type": "string", "required": False},
        "startToCloseTimeout": {"type": "string", "required": False},
        "taskList": {"type": "TaskList", "required": False},
        "taskPriority": {"type": "string", "required": False},
    },
    "scheduleLambdaFunctionDecisionAttributes": {
        "id": {"type": "string", "required": True},
        "input": {"type": "string", "required": False},
        "name": {"type": "string", "required": True},
        "startToCloseTimeout": {"type": "string", "required": False},
    },
    "signalExternalWorkflowExecutionDecisionAttributes": {
        "control": {"type": "string", "required": False},
        "input": {"type": "string", "required": False},
        "runId": {"type": "string", "required": False},
        "signalName": {"type": "string", "required": True},
        "workflowId": {"type": "string", "required": True},
    },
    "startChildWorkflowExecutionDecisionAttributes": {
        "childPolicy": {"type": "string", "required": False},
        "control": {"type": "string", "required": False},
        "executionStartToCloseTimeout": {"type": "string", "required": False},
        "input": {"type": "string", "required": False},
        "lambdaRole": {"type": "string", "required": False},
        "tagList": {"type": "string", "array": True, "required": False},
        "taskList": {"type": "TaskList", "required": False},
        "taskPriority": {"type": "string", "required": False},
        "taskStartToCloseTimeout": {"type": "string", "required": False},
        "workflowId": {"type": "string", "required": True},
        "workflowType": {"type": "WorkflowType", "required": True},
    },
    "startTimerDecisionAttributes": {
        "control": {"type": "string", "required": False},
        "startToFireTimeout": {"type": "string", "required": True},
        "timerId": {"type": "string", "required": True},
    },
}
