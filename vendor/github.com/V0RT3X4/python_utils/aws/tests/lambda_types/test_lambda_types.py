import unittest
from vortexa_utils.aws.lambdr.types import LambdaDict, LambdaContext
from .message_eg import handler as handler_message, get_message
from .repeat_eg import handler as handler_repeat, get_output


class TestMessageFunction(unittest.TestCase):

    def setUp(self):
        self.context = LambdaContext()

    def test_handler(self) -> None:
        event: LambdaDict = {
            "first_name": "Alex",
            "last_name": "Casalboni",
        }
        result = handler_message(event, self.context)
        self.assertIn('message', result)

    def test_handler_empty(self) -> None:
        event: LambdaDict = {}
        result = handler_message(event, self.context)
        self.assertIn('message', result)

    def test_message_default(self) -> None:
        msg = get_message()
        self.assertIsInstance(msg, str)
        self.assertIn('Hello', msg)
        self.assertIn('John', msg)
        self.assertIn('Smith', msg)
        self.assertTrue(msg.endswith('!'))

    def test_message_firstname(self) -> None:
        msg = get_message(first_name='Charlie')
        self.assertIsInstance(msg, str)
        self.assertIn('Hello', msg)
        self.assertIn('Charlie', msg)
        self.assertIn('Smith', msg)
        self.assertTrue(msg.endswith('!'))

    def test_message_lastname(self) -> None:
        msg = get_message(last_name='Brown')
        self.assertIsInstance(msg, str)
        self.assertIn('Hello', msg)
        self.assertIn('John', msg)
        self.assertIn('Brown', msg)
        self.assertTrue(msg.endswith('!'))

    def test_message(self) -> None:
        msg = get_message(first_name='Charlie', last_name='Brown')
        self.assertIsInstance(msg, str)
        self.assertIn('Hello', msg)
        self.assertIn('Charlie', msg)
        self.assertIn('Brown', msg)
        self.assertTrue(msg.endswith('!'))


class TestRepeatFunction(unittest.TestCase):

    def setUp(self):
        self.context = LambdaContext()

    def test_handler(self) -> None:
        event: LambdaDict = {
            "input": "NaN",
        }
        result = handler_repeat(event, self.context)
        self.assertIn('output', result)
        self.assertEqual(30, len(result['output']))

    def test_handler_empty(self) -> None:
        event: LambdaDict = {}
        with self.assertRaises(KeyError):
            handler_repeat(event, self.context)

    def test_repeat_empty_string(self) -> None:
        output = get_output('', 100)
        self.assertIsInstance(output, str)
        self.assertEqual(0, len(output))

    def test_repeat_zero(self) -> None:
        output = get_output('hello', 0)
        self.assertIsInstance(output, str)
        self.assertEqual(0, len(output))

    def test_repeat(self) -> None:
        output = get_output('hello', 10)
        self.assertIsInstance(output, str)
        self.assertEqual(50, len(output))
