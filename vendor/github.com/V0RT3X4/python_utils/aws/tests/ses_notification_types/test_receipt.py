from json import loads
from vortexa_utils.aws.ses.notification.types import Receipt


receipt_json = loads("""
{
"timestamp": "2015-09-11T20:32:33.936Z",
"processingTimeMillis": 222,
"recipients": [
	"recipient@example.com"
],
"spamVerdict": {
	"status": "PASS"
},
"virusVerdict": {
	"status": "PASS"
},
"spfVerdict": {
	"status": "PASS"
},
"dkimVerdict": {
	"status": "PASS"
},
"action": {
	"type": "SNS",
	"topicArn": "arn:aws:sns:us-east-1:012345678912:example-topic"
}
}
""")


def test_receipt():
    receipt = Receipt(**receipt_json)
    receipt.dkimVerdict.status == "PASS"
