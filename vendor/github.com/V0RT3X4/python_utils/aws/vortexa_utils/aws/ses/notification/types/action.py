from dataclasses import dataclass


@dataclass
class Action:
    """Action Object.

    Attributes
    ----------
    type : str
        action that was executed. [S3, SNS, Bounce, Lambda, Stop, WorkMail].
    topicArn : str
        Amazon Resource Name (ARN) of the SNS topic of the notification.
    bucketName : str
        S3 bucket to which the message was published.
        *Present only for the S3 action type.*
    objectKey : str
        name that uniquely identifies the email in the Amazon S3 bucket.
        This is the same as the messageId in the mail Object.
        *Present only for the S3 action type.*
    smtpReplyCode : str
         SMTP reply code, as defined by RFC 5321.
         *Present only for the bounce action type.*
    statusCode : str
        SMTP enhanced status code, as defined by RFC 3463.
        *Present only for the bounce action type.*
    message : str
        human-readable text to include in the bounce message.
        *Present only for the bounce action type.*
    sender : str
        The email address of the sender of the email that bounced.
        This is the address from which the bounce message was sent.
        *Present only for the bounce action type.*
    functionArn : str
        ARN of the Lambda function that was triggered.
        *Present only for the Lambda action type.*
    invocationType : str
        invocation type of the Lambda function. [RequestResponse, Event]
        *Present only for the Lambda action type.*
    organizationArn : str
         ARN of the Amazon WorkMail organization.
         *Present only for the WorkMail action type.*

    _see <https://docs.aws.amazon.com/ses/latest/DeveloperGuide/receiving-email-notifications-contents.html#receiving-email-notifications-contents-action-object>
    """
    type: str
    topicArn: str = None
    bucketName: str = None
    objectKey: str = None
    smtpReplyCode: str = None
    statusCode: str = None
    message: str = None
    sender: str = None
    functionArn: str = None
    invocationType: str = None
    organizationArn: str = None
