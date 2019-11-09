""" Example #2 """
import os
from vortexa_utils.aws.lambdr.types import LambdaDict, LambdaContext

N: int = int(os.environ.get('N') or 10)
STAGE: str = os.environ.get('stage') or 'dev'


def handler(event: LambdaDict, context: LambdaContext) -> LambdaDict:
    print('Received event {} for stage {}'.format(event, STAGE))
    input: str = event['input']  # required
    return {
        'output': get_output(input, N),
    }


def get_output(input: str, num: int):
    """ Return the input string repeated N times. """
    return input * num
