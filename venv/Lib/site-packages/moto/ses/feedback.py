"""
SES Feedback messages
Extracted from https://docs.aws.amazon.com/ses/latest/DeveloperGuide/notification-contents.html
"""

COMMON_MAIL = {
    "notificationType": "Bounce, Complaint, or Delivery.",
    "mail": {
        "timestamp": "2018-10-08T14:05:45 +0000",
        "messageId": "000001378603177f-7a5433e7-8edb-42ae-af10-f0181f34d6ee-000000",
        "source": "sender@example.com",
        "sourceArn": "arn:aws:ses:us-west-2:888888888888:identity/example.com",
        "sourceIp": "127.0.3.0",
        "sendingAccountId": None,
        "destination": ["recipient@example.com"],
        "headersTruncated": False,
        "headers": [
            {"name": "From", "value": '"Sender Name" <sender@example.com>'},
            {"name": "To", "value": '"Recipient Name" <recipient@example.com>'},
        ],
        "commonHeaders": {
            "from": ["Sender Name <sender@example.com>"],
            "date": "Mon, 08 Oct 2018 14:05:45 +0000",
            "to": ["Recipient Name <recipient@example.com>"],
            "messageId": " custom-message-ID",
            "subject": "Message sent using Amazon SES",
        },
    },
}
BOUNCE = {
    "bounceType": "Permanent",
    "bounceSubType": "General",
    "bouncedRecipients": [
        {
            "status": "5.0.0",
            "action": "failed",
            "diagnosticCode": "smtp; 550 user unknown",
            "emailAddress": "recipient1@example.com",
        },
        {
            "status": "4.0.0",
            "action": "delayed",
            "emailAddress": "recipient2@example.com",
        },
    ],
    "reportingMTA": "example.com",
    "timestamp": "2012-05-25T14:59:38.605Z",
    "feedbackId": "000001378603176d-5a4b5ad9-6f30-4198-a8c3-b1eb0c270a1d-000000",
    "remoteMtaIp": "127.0.2.0",
}
COMPLAINT = {
    "userAgent": "AnyCompany Feedback Loop (V0.01)",
    "complainedRecipients": [{"emailAddress": "recipient1@example.com"}],
    "complaintFeedbackType": "abuse",
    "arrivalDate": "2009-12-03T04:24:21.000-05:00",
    "timestamp": "2012-05-25T14:59:38.623Z",
    "feedbackId": "000001378603177f-18c07c78-fa81-4a58-9dd1-fedc3cb8f49a-000000",
}
DELIVERY = {
    "timestamp": "2014-05-28T22:41:01.184Z",
    "processingTimeMillis": 546,
    "recipients": ["success@simulator.amazonses.com"],
    "smtpResponse": "250 ok:  Message 64111812 accepted",
    "reportingMTA": "a8-70.smtp-out.amazonses.com",
    "remoteMtaIp": "127.0.2.0",
}
