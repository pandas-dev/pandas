import pytest

from pandas.io.formats.console import detect_console_encoding


class MockEncoding(object):  # TODO(py27): replace with mock
    """
    Used to add a side effect when accessing the 'encoding' property. If the
    side effect is a str in nature, the value will be returned. Otherwise, the
    side effect should be an exception that will be raised.
    """
    def __init__(self, encoding):
        super(MockEncoding, self).__init__()
        self.val = encoding

    @property
    def encoding(self):
        return raise_or_return(self.val)


def raise_or_return(val):
    if isinstance(val, str):
        return val
    else:
        raise val


@pytest.mark.parametrize('empty,filled', [
    ['stdin', 'stdout'],
    ['stdout', 'stdin']
])
def test_detect_console_encoding_from_stdout_stdin(monkeypatch, empty, filled):
    # Ensures that when sys.stdout.encoding or sys.stdin.encoding is used when
    # they have values filled.
    # GH 21552
    with monkeypatch.context() as context:
        context.setattr('sys.{}'.format(empty), MockEncoding(''))
        context.setattr('sys.{}'.format(filled), MockEncoding(filled))
        assert detect_console_encoding() == filled


@pytest.mark.parametrize('encoding', [
    MockEncoding(AttributeError),
    MockEncoding(IOError),
    MockEncoding('ascii')
])
def test_detect_console_encoding_fallback_to_locale(monkeypatch, encoding):
    # GH 21552
    with monkeypatch.context() as context:
        context.setattr('locale.getpreferredencoding', lambda: 'foo')
        context.setattr('sys.stdout', encoding)
        assert detect_console_encoding() == 'foo'


@pytest.mark.parametrize('std,locale', [
    [MockEncoding('ascii'), lambda: 'ascii'],
    [MockEncoding('ascii'), lambda: raise_or_return(Exception)],
    [MockEncoding(AttributeError), lambda: 'ascii'],
    [MockEncoding(AttributeError), lambda: raise_or_return(Exception)],
    [MockEncoding(IOError), lambda: 'ascii'],
    [MockEncoding(IOError), lambda: raise_or_return(Exception)]
])
def test_detect_console_encoding_fallback_to_default(monkeypatch, std, locale):
    # When both the stdout/stdin encoding and locale preferred encoding checks
    # fail (or return 'ascii', we should default to the sys default encoding.
    # GH 21552
    with monkeypatch.context() as context:
        context.setattr('locale.getpreferredencoding', locale)
        context.setattr('sys.stdout', std)
        context.setattr('sys.getdefaultencoding', lambda: 'sysDefaultEncoding')
        assert detect_console_encoding() == 'sysDefaultEncoding'
