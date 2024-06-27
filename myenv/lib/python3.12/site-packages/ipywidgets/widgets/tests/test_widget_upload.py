# Copyright (c) Jupyter Development Team.
# Distributed under the terms of the Modified BSD License.

import datetime as dt
from unittest import TestCase
from unittest.mock import MagicMock

from traitlets import TraitError

from ipywidgets import FileUpload


FILE_UPLOAD_FRONTEND_CONTENT = {
    'name': 'file-name.txt',
    'type': 'text/plain',
    'size': 20760,
    'last_modified': 1578578296434,
    'content': memoryview(b'file content'),
}


class TestFileUpload(TestCase):

    def test_construction(self):
        uploader = FileUpload()
        # Default
        assert uploader.accept == ''
        assert not uploader.multiple
        assert not uploader.disabled

    def test_construction_with_params(self):
        uploader = FileUpload(
            accept='.txt', multiple=True, disabled=True)
        assert uploader.accept == '.txt'
        assert uploader.multiple
        assert uploader.disabled

    def test_empty_initial_value(self):
        uploader = FileUpload()
        assert uploader.value == ()

    def test_receive_single_file(self):
        uploader = FileUpload()
        message = {'value': [FILE_UPLOAD_FRONTEND_CONTENT]}
        uploader.set_state(message)
        assert len(uploader.value) == 1
        (uploaded_file,) = uploader.value
        assert uploaded_file.name == 'file-name.txt'
        assert uploaded_file.type == 'text/plain'
        assert uploaded_file.size == 20760
        assert uploaded_file.content.tobytes() == b'file content'
        assert (
            uploaded_file.last_modified ==
            dt.datetime(2020, 1, 9, 13, 58, 16, 434000, tzinfo=dt.timezone.utc)
        )

    def test_receive_multiple_files(self):
        uploader = FileUpload(multiple=True)
        message = {
            'value': [
                FILE_UPLOAD_FRONTEND_CONTENT,
                {**FILE_UPLOAD_FRONTEND_CONTENT, **{'name': 'other-file-name.txt'}}
            ]
        }
        uploader.set_state(message)
        assert len(uploader.value) == 2
        assert uploader.value[0].name == 'file-name.txt'
        assert uploader.value[1].name == 'other-file-name.txt'

    def test_serialization_deserialization_integrity(self):
        # The value traitlet needs to remain unchanged following
        # a serialization / deserialization roundtrip, otherwise
        # the kernel dispatches it back to the frontend following
        # a state change, because it doesn't recognize that the
        # property_lock entry is the same as the new value.
        from ipykernel.comm import Comm
        uploader = FileUpload()
        mock_comm = MagicMock(spec=Comm)
        mock_comm.send = MagicMock()
        mock_comm.kernel = 'does not matter'
        uploader.comm = mock_comm
        message = {'value': [FILE_UPLOAD_FRONTEND_CONTENT]}
        uploader.set_state(message)

        # Check that no message is sent back to the frontend
        # as a result of setting the state.
        mock_comm.send.assert_not_called()

    def test_resetting_value(self):
        # Simulate an upload, then resetting the value from the
        # kernel.
        uploader = FileUpload()
        message = {'value': [FILE_UPLOAD_FRONTEND_CONTENT]}
        uploader.set_state(message)

        uploader.value = []  # reset value to an empty file list

        assert uploader.get_state(key='value') == {'value': []}

    def test_setting_non_empty_value(self):
        # Simulate user setting a value for the upload from the kernel.
        uploader = FileUpload()
        content = memoryview(b'some content')
        uploader.value = [{
            'name': 'some-name.txt',
            'type': 'text/plain',
            'size': 561,
            'last_modified': dt.datetime(2020, 1, 9, 13, 58, 16, 434000, tzinfo=dt.timezone.utc),
            'content': content
        }]
        state = uploader.get_state(key='value')
        assert len(state['value']) == 1
        [entry] = state['value']
        assert entry['name'] == 'some-name.txt'
        assert entry['type'] == 'text/plain'
        assert entry['size'] == 561
        assert entry['last_modified'] == 1578578296434
        assert entry['content'] == content
