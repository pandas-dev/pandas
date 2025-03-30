# Copyright(c) Jupyter Development Team.
# Distributed under the terms of the Modified BSD License.

"""FileUpload class.

Represents a file upload button.
"""
import datetime as dt

from traitlets import (
    observe, default, Unicode, Dict, Int, Bool, Bytes, CaselessStrEnum
)

from .widget_description import DescriptionWidget
from .valuewidget import ValueWidget
from .widget_core import CoreWidget
from .widget_button import ButtonStyle
from .widget import register, widget_serialization
from .trait_types import InstanceDict, TypedTuple
from traitlets import Bunch


def _deserialize_single_file(js):
    uploaded_file = Bunch()
    for attribute in ['name', 'type', 'size', 'content']:
        uploaded_file[attribute] = js[attribute]
    uploaded_file['last_modified'] = dt.datetime.fromtimestamp(
        js['last_modified'] / 1000,
        tz=dt.timezone.utc
    )
    return uploaded_file


def _deserialize_value(js, _):
    return [_deserialize_single_file(entry) for entry in js]


def _serialize_single_file(uploaded_file):
    js = {}
    for attribute in ['name', 'type', 'size', 'content']:
        js[attribute] = uploaded_file[attribute]
    js['last_modified'] = int(uploaded_file['last_modified'].timestamp() * 1000)
    return js


def _serialize_value(value, _):
    return [_serialize_single_file(entry) for entry in value]


_value_serialization = {
    'from_json': _deserialize_value,
    'to_json': _serialize_value
}


@register
class FileUpload(DescriptionWidget, ValueWidget, CoreWidget):
    """File upload widget

    This creates a file upload input that allows the user to select
    one or more files to upload. The file metadata and content
    can be retrieved in the kernel.

    Examples
    --------

    >>> import ipywidgets as widgets
    >>> uploader = widgets.FileUpload()

    # After displaying `uploader` and uploading a file:

    >>> uploader.value
    [
      {
        'name': 'example.txt',
        'type': 'text/plain',
        'size': 36,
        'last_modified': datetime.datetime(2020, 1, 9, 15, 58, 43, 321000, tzinfo=datetime.timezone.utc),
        'content': <memory at 0x10c1b37c8>
      }
    ]
    >>> uploader.value[0].content.tobytes()
    b'This is the content of example.txt.\n'

    Parameters
    ----------

    accept: str, optional
        Which file types to accept, e.g. '.doc,.docx'. For a full
        description of how to specify this, see
        https://developer.mozilla.org/en-US/docs/Web/HTML/Element/input/file#attr-accept
        Defaults to accepting all file types.

    multiple: bool, optional
        Whether to accept multiple files at the same time. Defaults to False.

    disabled: bool, optional
        Whether user interaction is enabled.

    icon: str, optional
        The icon to use for the button displayed on the screen.
        Can be any Font-awesome icon without the fa- prefix.
        Defaults to 'upload'. If missing, no icon is shown.

    description: str, optional
        The text to show on the label. Defaults to 'Upload'.

    button_style: str, optional
        One of 'primary', 'success', 'info', 'warning', 'danger' or ''.

    style: widgets.widget_button.ButtonStyle, optional
        Style configuration for the button.

    value: Tuple[Dict], optional
        The value of the last uploaded file or set of files. See the
        documentation for details of how to use this to retrieve file
        content and metadata:
        https://ipywidgets.readthedocs.io/en/stable/examples/Widget%20List.html#File-Upload

    error: str, optional
        Whether the last upload triggered an error.
    """
    _model_name = Unicode('FileUploadModel').tag(sync=True)
    _view_name = Unicode('FileUploadView').tag(sync=True)

    accept = Unicode(help='File types to accept, empty string for all').tag(sync=True)
    multiple = Bool(help='If True, allow for multiple files upload').tag(sync=True)
    disabled = Bool(help='Enable or disable button').tag(sync=True)
    icon = Unicode('upload', help="Font-awesome icon name, without the 'fa-' prefix.").tag(sync=True)
    button_style = CaselessStrEnum(
        values=['primary', 'success', 'info', 'warning', 'danger', ''], default_value='',
        help='Use a predefined styling for the button.').tag(sync=True)
    style = InstanceDict(ButtonStyle).tag(sync=True, **widget_serialization)
    error = Unicode(help='Error message').tag(sync=True)
    value = TypedTuple(Dict(), help='The file upload value').tag(
        sync=True, echo_update=False, **_value_serialization)

    @default('description')
    def _default_description(self):
        return 'Upload'
