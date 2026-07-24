import sys
from unittest import TestCase
from contextlib import contextmanager

from IPython.display import Markdown, Image
from ipywidgets import widget_output


class TestOutputWidget(TestCase):

    @contextmanager
    def _mocked_ipython(self, get_ipython, clear_output):
        """ Context manager that monkeypatches get_ipython and clear_output """
        original_clear_output = widget_output.clear_output
        original_get_ipython = widget_output.get_ipython
        widget_output.get_ipython = get_ipython
        widget_output.clear_output = clear_output
        try:
            yield
        finally:
            widget_output.clear_output = original_clear_output
            widget_output.get_ipython = original_get_ipython

    def _mock_get_ipython(self, msg_id):
        """ Returns a mock IPython application with a mocked kernel """
        kernel = type(
            'mock_kernel',
            (object, ),
            {'_parent_header': {'header': {'msg_id': msg_id}}}
        )

        # Specifically override this so the traceback
        # is still printed to screen
        def showtraceback(self_, exc_tuple, *args, **kwargs):
            etype, evalue, tb = exc_tuple
            raise etype(evalue)

        ipython = type(
            'mock_ipython',
            (object, ),
            {'kernel': kernel, 'showtraceback': showtraceback}
        )
        return ipython

    def _mock_clear_output(self):
        """ Mock function that records calls to it """
        calls = []

        def clear_output(*args, **kwargs):
            calls.append((args, kwargs))
        clear_output.calls = calls

        return clear_output

    def test_set_msg_id_when_capturing(self):
        msg_id = 'msg-id'
        get_ipython = self._mock_get_ipython(msg_id)
        clear_output = self._mock_clear_output()

        with self._mocked_ipython(get_ipython, clear_output):
            widget = widget_output.Output()
            assert widget.msg_id == ''
            with widget:
                assert widget.msg_id == msg_id
            assert widget.msg_id == ''

    def test_clear_output(self):
        msg_id = 'msg-id'
        get_ipython = self._mock_get_ipython(msg_id)
        clear_output = self._mock_clear_output()

        with self._mocked_ipython(get_ipython, clear_output):
            widget = widget_output.Output()
            widget.clear_output(wait=True)

        assert len(clear_output.calls) == 1
        assert clear_output.calls[0] == ((), {'wait': True})

    def test_capture_decorator(self):
        msg_id = 'msg-id'
        get_ipython = self._mock_get_ipython(msg_id)
        clear_output = self._mock_clear_output()
        expected_argument = 'arg'
        expected_keyword_argument = True
        captee_calls = []

        with self._mocked_ipython(get_ipython, clear_output):
            widget = widget_output.Output()
            assert widget.msg_id == ''

            @widget.capture()
            def captee(*args, **kwargs):
                # Check that we are capturing output
                assert widget.msg_id == msg_id

                # Check that arguments are passed correctly
                captee_calls.append((args, kwargs))

            captee(
                expected_argument, keyword_argument=expected_keyword_argument)
            assert widget.msg_id == ''
            captee()

        assert len(captee_calls) == 2
        assert captee_calls[0] == (
            (expected_argument, ),
            {'keyword_argument': expected_keyword_argument}
        )
        assert captee_calls[1] == ((), {})

    def test_capture_decorator_clear_output(self):
        msg_id = 'msg-id'
        get_ipython = self._mock_get_ipython(msg_id)
        clear_output = self._mock_clear_output()

        with self._mocked_ipython(get_ipython, clear_output):
            widget = widget_output.Output()

            @widget.capture(clear_output=True, wait=True)
            def captee(*args, **kwargs):
                # Check that we are capturing output
                assert widget.msg_id == msg_id

            captee()
            captee()

        assert len(clear_output.calls) == 2
        assert clear_output.calls[0] == clear_output.calls[1] == \
            ((), {'wait': True})

    def test_capture_decorator_no_clear_output(self):
        msg_id = 'msg-id'
        get_ipython = self._mock_get_ipython(msg_id)
        clear_output = self._mock_clear_output()

        with self._mocked_ipython(get_ipython, clear_output):
            widget = widget_output.Output()

            @widget.capture(clear_output=False)
            def captee(*args, **kwargs):
                # Check that we are capturing output
                assert widget.msg_id == msg_id

            captee()
            captee()

        assert len(clear_output.calls) == 0


def _make_stream_output(text, name):
    return {
        'output_type': 'stream',
        'name': name,
        'text': text
    }


def test_append_stdout():
    widget = widget_output.Output()

    # Try appending a message to stdout.
    widget.append_stdout("snakes!")
    expected = (_make_stream_output("snakes!", "stdout"),)
    assert widget.outputs == expected, repr(widget.outputs)

    # Try appending a second message.
    widget.append_stdout("more snakes!")
    expected += (_make_stream_output("more snakes!", "stdout"),)
    assert widget.outputs == expected, repr(widget.outputs)


def test_append_stderr():
    widget = widget_output.Output()

    # Try appending a message to stderr.
    widget.append_stderr("snakes!")
    expected = (_make_stream_output("snakes!", "stderr"),)
    assert widget.outputs == expected, repr(widget.outputs)

    # Try appending a second message.
    widget.append_stderr("more snakes!")
    expected += (_make_stream_output("more snakes!", "stderr"),)
    assert widget.outputs == expected, repr(widget.outputs)


def test_append_display_data():
    widget = widget_output.Output()

    # Try appending a Markdown object.
    widget.append_display_data(Markdown("# snakes!"))
    expected = (
        {
            'output_type': 'display_data',
            'data': {
                'text/plain': '<IPython.core.display.Markdown object>',
                'text/markdown': '# snakes!'
            },
            'metadata': {}
        },
    )
    assert widget.outputs == expected, repr(widget.outputs)

    # Now try appending an Image.
    image_data = b"foobar"

    widget.append_display_data(Image(image_data, width=123, height=456))
    # Old ipykernel/IPython
    expected1 = expected + (
        {
            'output_type': 'display_data',
            'data': {
                'image/png': 'Zm9vYmFy\n',
                'text/plain': '<IPython.core.display.Image object>'
            },
            'metadata': {
                'image/png': {
                    'width': 123,
                    'height': 456
                }
            }
        },
    )
    # Latest ipykernel/IPython
    expected2 = expected + (
        {
            'output_type': 'display_data',
            'data': {
                'image/png': 'Zm9vYmFy',
                'text/plain': '<IPython.core.display.Image object>'
            },
            'metadata': {
                'image/png': {
                    'width': 123,
                    'height': 456
                }
            }
        },
    )
    assert widget.outputs == expected1 or widget.outputs == expected2
