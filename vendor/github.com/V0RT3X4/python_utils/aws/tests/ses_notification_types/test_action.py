from vortexa_utils.aws.ses.notification.types import Action
from json import loads


action_json_sns = """
{
	"type": "SNS",
	"topicArn": "arn:aws:sns:us-east-1:012345678912:example-topic"
}
"""


def test_sns_action():
    action = Action(**loads(action_json_sns))
    assert action.type == "SNS"
    assert action.topicArn == "arn:aws:sns:us-east-1:012345678912:example-topic"
