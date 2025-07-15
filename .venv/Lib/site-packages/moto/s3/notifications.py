import copy
import json
from datetime import datetime
from enum import Enum
from typing import TYPE_CHECKING, Any, Dict, List

from moto.core.utils import unix_time

if TYPE_CHECKING:
    from moto.s3.models import FakeBucket

_EVENT_TIME_FORMAT = "%Y-%m-%dT%H:%M:%S.%f"


class S3NotificationEvent(str, Enum):
    REDUCED_REDUNDANCY_LOST_OBJECT_EVENT = "s3:ReducedRedundancyLostObject"
    OBJECT_CREATED_EVENT = "s3:ObjectCreated:*"
    OBJECT_CREATED_PUT_EVENT = "s3:ObjectCreated:Put"
    OBJECT_CREATED_POST_EVENT = "s3:ObjectCreated:Post"
    OBJECT_CREATED_COPY_EVENT = "s3:ObjectCreated:Copy"
    OBJECT_CREATED_COMPLETE_MULTIPART_UPLOAD_EVENT = (
        "s3:ObjectCreated:CompleteMultipartUpload"
    )
    OBJECT_REMOVED_EVENT = "s3:ObjectRemoved:*"
    OBJECT_REMOVED_DELETE_EVENT = "s3:ObjectRemoved:Delete"
    OBJECT_REMOVED_DELETE_MARKER_CREATED_EVENT = "s3:ObjectRemoved:DeleteMarkerCreated"
    OBJECT_RESTORE_EVENT = "s3:ObjectRestore:*"
    OBJECT_RESTORE_POST_EVENT = "s3:ObjectRestore:Post"
    OBJECT_RESTORE_COMPLETED_EVENT = "s3:ObjectRestore:Completed"
    REPLICATION_EVENT = "s3:Replication:*"
    REPLICATION_OPERATION_FAILED_REPLICATION_EVENT = (
        "s3:Replication:OperationFailedReplication"
    )
    REPLICATION_OPERATION_NOT_TRACKED_EVENT = "s3:Replication:OperationNotTracked"
    REPLICATION_OPERATION_MISSED_THRESHOLD_EVENT = (
        "s3:Replication:OperationMissedThreshold"
    )
    REPLICATION_OPERATION_REPLICATED_AFTER_THRESHOLD_EVENT = (
        "s3:Replication:OperationReplicatedAfterThreshold"
    )
    OBJECT_RESTORE_DELETE_EVENT = "s3:ObjectRestore:Delete"
    LIFECYCLE_TRANSITION_EVENT = "s3:LifecycleTransition"
    INTELLIGENT_TIERING_EVENT = "s3:IntelligentTiering"
    OBJECT_ACL_UPDATE_EVENT = "s3:ObjectAcl:Put"
    LIFECYCLE_EXPIRATION_EVENT = "s3:LifecycleExpiration:*"
    LIFECYCLEEXPIRATION_DELETE_EVENT = "s3:LifecycleExpiration:Delete"
    LIFECYCLE_EXPIRATION_DELETE_MARKER_CREATED_EVENT = (
        "s3:LifecycleExpiration:DeleteMarkerCreated"
    )
    OBJECT_TAGGING_EVENT = "s3:ObjectTagging:*"
    OBJECT_TAGGING_PUT_EVENT = "s3:ObjectTagging:Put"
    OBJECT_TAGGING_DELETE_EVENT = "s3:ObjectTagging:Delete"

    @classmethod
    def events(self) -> List[str]:
        return sorted([item.value for item in S3NotificationEvent])

    @classmethod
    def is_event_valid(self, event_name: str) -> bool:
        # Ex) s3:ObjectCreated:Put
        if event_name in self.events():
            return True
        # Ex) event name without `s3:` like ObjectCreated:Put
        if event_name in [e[:3] for e in self.events()]:
            return True
        return False


def _get_s3_event(
    event_name: str, bucket: "FakeBucket", key: Any, notification_id: str
) -> Dict[str, List[Dict[str, Any]]]:
    etag = key.etag.replace('"', "")
    # s3:ObjectCreated:Put --> ObjectCreated:Put
    event_name = event_name[3:]
    event_time = datetime.now().strftime(_EVENT_TIME_FORMAT)
    return {
        "Records": [
            {
                "eventVersion": "2.1",
                "eventSource": "aws:s3",
                "awsRegion": bucket.region_name,
                "eventTime": event_time,
                "eventName": event_name,
                "s3": {
                    "s3SchemaVersion": "1.0",
                    "configurationId": notification_id,
                    "bucket": {
                        "name": bucket.name,
                        "arn": bucket.arn,
                    },
                    "object": {"key": key.name, "size": key.size, "eTag": etag},
                },
            }
        ]
    }


def _get_region_from_arn(arn: str) -> str:
    return arn.split(":")[3]


def send_event(
    account_id: str, event_name: S3NotificationEvent, bucket: Any, key: Any
) -> None:
    if bucket.notification_configuration is None:
        return

    for notification in bucket.notification_configuration.cloud_function:
        if notification.matches(event_name, key.name):
            event_body = _get_s3_event(event_name, bucket, key, notification.id)
            region_name = _get_region_from_arn(notification.arn)

            _invoke_awslambda(account_id, event_body, notification.arn, region_name)

    for notification in bucket.notification_configuration.queue:
        if notification.matches(event_name, key.name):
            event_body = _get_s3_event(event_name, bucket, key, notification.id)
            region_name = _get_region_from_arn(notification.arn)
            queue_name = notification.arn.split(":")[-1]

            _send_sqs_message(account_id, event_body, queue_name, region_name)

    for notification in bucket.notification_configuration.topic:
        if notification.matches(event_name, key.name):
            event_body = _get_s3_event(event_name, bucket, key, notification.id)
            region_name = _get_region_from_arn(notification.arn)
            topic_arn = notification.arn

            _send_sns_message(account_id, event_body, topic_arn, region_name)

    if bucket.notification_configuration.event_bridge is not None:
        _send_event_bridge_message(account_id, bucket, event_name, key)


def _send_sqs_message(
    account_id: str, event_body: Any, queue_name: str, region_name: str
) -> None:
    try:
        from moto.sqs.models import sqs_backends

        sqs_backend = sqs_backends[account_id][region_name]
        sqs_backend.send_message(
            queue_name=queue_name, message_body=json.dumps(event_body)
        )
    except:  # noqa
        # This is an async action in AWS.
        # Even if this part fails, the calling function should pass, so catch all errors
        # Possible exceptions that could be thrown:
        # - Queue does not exist
        pass


def _send_sns_message(
    account_id: str, event_body: Any, topic_arn: str, region_name: str
) -> None:
    try:
        from moto.sns.models import sns_backends

        sns_backend = sns_backends[account_id][region_name]
        sns_backend.publish(arn=topic_arn, message=json.dumps(event_body))
    except:  # noqa
        # This is an async action in AWS.
        # Even if this part fails, the calling function should pass, so catch all errors
        # Possible exceptions that could be thrown:
        # - Topic does not exist
        pass


def _send_event_bridge_message(
    account_id: str,
    bucket: "FakeBucket",
    event_name: str,
    key: Any,
) -> None:
    try:
        from moto.events.models import events_backends
        from moto.events.utils import _BASE_EVENT_MESSAGE

        event = copy.deepcopy(_BASE_EVENT_MESSAGE)
        event["detail-type"] = _detail_type(event_name)
        event["source"] = "aws.s3"
        event["account"] = account_id
        event["time"] = unix_time()
        event["region"] = bucket.region_name
        event["resources"] = [bucket.arn]
        event["detail"] = {
            "version": "0",
            "bucket": {"name": bucket.name},
            "object": {
                "key": key.name,
                "size": key.size,
                "eTag": key.etag.replace('"', ""),
                "version-id": key.version_id,
                "sequencer": "617f08299329d189",
            },
            "request-id": "N4N7GDK58NMKJ12R",
            "requester": "123456789012",
            "source-ip-address": "1.2.3.4",
            # ex) s3:ObjectCreated:Put -> ObjectCreated
            "reason": event_name.split(":")[1],
        }

        events_backend = events_backends[account_id][bucket.region_name]
        for event_bus in events_backend.event_buses.values():
            for rule in event_bus.rules.values():
                rule.send_to_targets(event, transform_input=False)

    except:  # noqa
        # This is an async action in AWS.
        # Even if this part fails, the calling function should pass, so catch all errors
        # Possible exceptions that could be thrown:
        # - EventBridge does not exist
        pass


def _detail_type(event_name: str) -> str:
    """Detail type field values for event messages of s3 EventBridge notification

    document: https://docs.aws.amazon.com/AmazonS3/latest/userguide/EventBridge.html
    """
    if event_name in [e for e in S3NotificationEvent.events() if "ObjectCreated" in e]:
        return "Object Created"
    elif event_name in [
        e
        for e in S3NotificationEvent.events()
        if "ObjectRemoved" in e or "LifecycleExpiration" in e
    ]:
        return "Object Deleted"
    elif event_name in [
        e for e in S3NotificationEvent.events() if "ObjectRestore" in e
    ]:
        if event_name == S3NotificationEvent.OBJECT_RESTORE_POST_EVENT:
            return "Object Restore Initiated"
        elif event_name == S3NotificationEvent.OBJECT_RESTORE_COMPLETED_EVENT:
            return "Object Restore Completed"
        else:
            # s3:ObjectRestore:Delete event
            return "Object Restore Expired"
    elif event_name in [
        e for e in S3NotificationEvent.events() if "LifecycleTransition" in e
    ]:
        return "Object Storage Class Changed"
    elif event_name in [
        e for e in S3NotificationEvent.events() if "IntelligentTiering" in e
    ]:
        return "Object Access Tier Changed"
    elif event_name in [e for e in S3NotificationEvent.events() if "ObjectAcl" in e]:
        return "Object ACL Updated"
    elif event_name in [e for e in S3NotificationEvent.events() if "ObjectTagging"]:
        if event_name == S3NotificationEvent.OBJECT_TAGGING_PUT_EVENT:
            return "Object Tags Added"
        else:
            # s3:ObjectTagging:Delete event
            return "Object Tags Deleted"
    else:
        raise ValueError(
            f"unsupported event `{event_name}` for s3 eventbridge notification (https://docs.aws.amazon.com/AmazonS3/latest/userguide/EventBridge.html)"
        )


def _invoke_awslambda(
    account_id: str, event_body: Any, fn_arn: str, region_name: str
) -> None:
    try:
        from moto.awslambda.utils import get_backend

        lambda_backend = get_backend(account_id, region_name)
        func = lambda_backend.get_function(fn_arn)
        func.invoke(json.dumps(event_body), dict(), dict())
    except:  # noqa
        # This is an async action in AWS.
        # Even if this part fails, the calling function should pass, so catch all errors
        # Possible exceptions that could be thrown:
        # - Function does not exist
        pass


def _get_test_event(bucket_name: str) -> Dict[str, Any]:
    event_time = datetime.now().strftime(_EVENT_TIME_FORMAT)
    return {
        "Service": "Amazon S3",
        "Event": "s3:TestEvent",
        "Time": event_time,
        "Bucket": bucket_name,
    }


def send_test_event(account_id: str, bucket: Any) -> None:
    arns = [n.arn for n in bucket.notification_configuration.queue]
    for arn in set(arns):
        region_name = _get_region_from_arn(arn)
        queue_name = arn.split(":")[-1]
        message_body = _get_test_event(bucket.name)
        _send_sqs_message(account_id, message_body, queue_name, region_name)

    arns = [n.arn for n in bucket.notification_configuration.topic]
    for arn in set(arns):
        region_name = _get_region_from_arn(arn)
        message_body = _get_test_event(bucket.name)
        _send_sns_message(account_id, message_body, arn, region_name)
