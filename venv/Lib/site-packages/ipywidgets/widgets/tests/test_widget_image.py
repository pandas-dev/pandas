# Copyright (c) Jupyter Development Team.
# Distributed under the terms of the Modified BSD License.

"""Test Image widget"""

import io
import os

from ipywidgets import Image

import hashlib

import pkgutil

import tempfile
from contextlib import contextmanager

# Data
@contextmanager
def get_logo_png():
    # Once the tests are not in the package, this context manager can be
    # replaced with the location of the actual file
    LOGO_DATA = pkgutil.get_data('ipywidgets.widgets.tests',
                                 'data/jupyter-logo-transparent.png')
    handle, fname = tempfile.mkstemp()
    os.close(handle)
    with open(fname, 'wb') as f:
        f.write(LOGO_DATA)

    yield fname

    os.remove(fname)

LOGO_PNG_DIGEST = '3ff9eafd7197083153e83339a72e7a335539bae189c33554c680e4382c98af02'


def test_empty_image():
    # Empty images shouldn't raise any errors
    Image()


def test_image_value():
    random_bytes = b'\x0ee\xca\x80\xcd\x9ak#\x7f\x07\x03\xa7'

    Image(value=random_bytes)


def test_image_format():
    # Test that these format names don't throw an error
    Image(format='png')

    Image(format='jpeg')

    Image(format='url')


def test_from_filename():
    with get_logo_png() as LOGO_PNG:
        img = Image.from_file(LOGO_PNG)

        assert_equal_hash(img.value, LOGO_PNG_DIGEST)


def test_set_from_filename():
    img = Image()
    with get_logo_png() as LOGO_PNG:
        img.set_value_from_file(LOGO_PNG)

        assert_equal_hash(img.value, LOGO_PNG_DIGEST)


def test_from_file():
    with get_logo_png() as LOGO_PNG:
        with open(LOGO_PNG, 'rb') as f:
            img = Image.from_file(f)
            assert_equal_hash(img.value, LOGO_PNG_DIGEST)


def test_set_value_from_file():
    img = Image()
    with get_logo_png() as LOGO_PNG:
        with open(LOGO_PNG, 'rb') as f:
            img.set_value_from_file(f)
            assert_equal_hash(img.value, LOGO_PNG_DIGEST)


def test_from_url_unicode():
    img = Image.from_url('https://jupyter.org/assets/main-logo.svg')
    assert img.value == b'https://jupyter.org/assets/main-logo.svg'


def test_from_url_bytes():
    img = Image.from_url(b'https://jupyter.org/assets/main-logo.svg')

    assert img.value == b'https://jupyter.org/assets/main-logo.svg'


def test_format_inference_filename():
    with tempfile.NamedTemporaryFile(suffix='.svg', delete=False) as f:
        name = f.name
        f.close()           # Allow tests to run on Windows
        img = Image.from_file(name)

    assert img.format == 'svg+xml'


def test_format_inference_file():
    with tempfile.NamedTemporaryFile(suffix='.gif', delete=False) as f:
        img = Image.from_file(f)

        assert img.format == 'gif'


def test_format_inference_stream():
    # There's no way to infer the format, so it should default to png
    fstream = io.BytesIO(b'')
    img = Image.from_file(fstream)

    assert img.format == 'png'


def test_serialize():
    fstream = io.BytesIO(b'123')
    img = Image.from_file(fstream)

    img_state = img.get_state()

    # for python27 it is a memoryview
    assert isinstance(img_state['value'], (bytes, memoryview))
    # make sure it is (for python 3), since that is what it will be once it comes off the wire
    img_state['value'] = memoryview(img_state['value'])

    # check that we can deserialize it and get back the original value
    img_copy = Image()
    img_copy.set_state(img_state)
    assert img.value == img_copy.value


def test_format_inference_overridable():
    with tempfile.NamedTemporaryFile(suffix='.svg', delete=False) as f:
        name = f.name
        f.close()           # Allow tests to run on Windows
        img = Image.from_file(name, format='gif')

    assert img.format == 'gif'


def test_value_repr_length():
    with get_logo_png() as LOGO_PNG:
        with open(LOGO_PNG, 'rb') as f:
            img = Image.from_file(f)
            assert len(img.__repr__()) < 140
            assert img.__repr__().endswith(")")
            assert img.__repr__()[-5:-2] == '...'


def test_value_repr_url():
    img = Image.from_url(b'https://jupyter.org/assets/main-logo.svg')

    assert 'https://jupyter.org/assets/main-logo.svg' in img.__repr__()


# Helper functions
def get_hash_hex(byte_str):
    m = hashlib.new('sha256')

    m.update(byte_str)

    return m.hexdigest()


def assert_equal_hash(byte_str, digest):
    assert get_hash_hex(byte_str) == digest
