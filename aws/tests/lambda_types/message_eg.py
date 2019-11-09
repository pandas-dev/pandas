""" Example #1 """
import os
from vortexa_utils.aws.lambdr.types import LambdaDict, LambdaContext

MSG_TEMPLATE: str = os.environ.get('MSG_TEMPLATE') or 'Hello {} {}!'
STAGE: str = os.environ.get('stage') or 'dev'


def handler(event: LambdaDict, context: LambdaContext) -> LambdaDict:
    print('Received event {} for stage {}'.format(event, STAGE))
    first_name: str = event.get('first_name')  # optional
    last_name: str = event.get('last_name')  # optional
    return {
        'message': get_message(first_name, last_name),
    }


def get_message(first_name: str = 'John', last_name: str = 'Smith'):
    return MSG_TEMPLATE.format(first_name, last_name)
