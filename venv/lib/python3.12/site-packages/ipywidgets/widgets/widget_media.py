# Copyright (c) Jupyter Development Team.
# Distributed under the terms of the Modified BSD License.

import mimetypes

from .widget_core import CoreWidget
from .domwidget import DOMWidget
from .valuewidget import ValueWidget
from .widget import register
from traitlets import Unicode, CUnicode, Bool
from .trait_types import CByteMemoryView


@register
class _Media(DOMWidget, ValueWidget, CoreWidget):
    """Base class for Image, Audio and Video widgets.

    The `value` of this widget accepts a byte string.  The byte string is the
    raw data that you want the browser to display.

    If you pass `"url"` to the `"format"` trait, `value` will be interpreted
    as a URL as bytes encoded in UTF-8.
    """

    # Define the custom state properties to sync with the front-end
    value = CByteMemoryView(help="The media data as a memory view of bytes.").tag(sync=True)

    @classmethod
    def _from_file(cls, tag, filename, **kwargs):
        """
        Create an :class:`Media` from a local file.

        Parameters
        ----------
        filename: str
            The location of a file to read into the value from disk.

        **kwargs:
            The keyword arguments for `Media`

        Returns an `Media` with the value set from the filename.
        """
        value = cls._load_file_value(filename)

        if 'format' not in kwargs:
            format = cls._guess_format(tag, filename)
            if format is not None:
                kwargs['format'] = format

        return cls(value=value, **kwargs)

    @classmethod
    def from_url(cls, url, **kwargs):
        """
        Create an :class:`Media` from a URL.

        :code:`Media.from_url(url)` is equivalent to:

        .. code-block: python

            med = Media(value=url, format='url')

        But both unicode and bytes arguments are allowed for ``url``.

        Parameters
        ----------
        url: [str, bytes]
            The location of a URL to load.
        """
        if isinstance(url, str):
            # If str, it needs to be encoded to bytes
            url = url.encode('utf-8')

        return cls(value=url, format='url', **kwargs)

    def set_value_from_file(self, filename):
        """
        Convenience method for reading a file into `value`.

        Parameters
        ----------
        filename: str
            The location of a file to read into value from disk.
        """
        value = self._load_file_value(filename)

        self.value = value

    @classmethod
    def _load_file_value(cls, filename):
        if getattr(filename, 'read', None) is not None:
            return filename.read()
        else:
            with open(filename, 'rb') as f:
                return f.read()

    @classmethod
    def _guess_format(cls, tag, filename):
        # file objects may have a .name parameter
        name = getattr(filename, 'name', None)
        name = name or filename

        try:
            mtype, _ = mimetypes.guess_type(name)
            if not mtype.startswith('{}/'.format(tag)):
                return None

            return mtype[len('{}/'.format(tag)):]
        except Exception:
            return None

    def _get_repr(self, cls):
        # Truncate the value in the repr, since it will
        # typically be very, very large.
        class_name = self.__class__.__name__

        # Return value first like a ValueWidget
        signature = []

        sig_value = 'value={!r}'.format(self.value[:40].tobytes())
        if self.value.nbytes > 40:
            sig_value = sig_value[:-1]+"..."+sig_value[-1]
        signature.append(sig_value)

        for key in super(cls, self)._repr_keys():
            if key == 'value':
                continue
            value = str(getattr(self, key))
            signature.append('{}={!r}'.format(key, value))
        signature = ', '.join(signature)
        return '{}({})'.format(class_name, signature)


@register
class Image(_Media):
    """Displays an image as a widget.

    The `value` of this widget accepts a byte string.  The byte string is the
    raw image data that you want the browser to display.  You can explicitly
    define the format of the byte string using the `format` trait (which
    defaults to "png").

    If you pass `"url"` to the `"format"` trait, `value` will be interpreted
    as a URL as bytes encoded in UTF-8.
    """
    _view_name = Unicode('ImageView').tag(sync=True)
    _model_name = Unicode('ImageModel').tag(sync=True)

    # Define the custom state properties to sync with the front-end
    format = Unicode('png', help="The format of the image.").tag(sync=True)
    width = CUnicode(help="Width of the image in pixels. Use layout.width "
                          "for styling the widget.").tag(sync=True)
    height = CUnicode(help="Height of the image in pixels. Use layout.height "
                           "for styling the widget.").tag(sync=True)

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

    @classmethod
    def from_file(cls, filename, **kwargs):
        return cls._from_file('image', filename, **kwargs)

    def __repr__(self):
        return self._get_repr(Image)


@register
class Video(_Media):
    """Displays a video as a widget.

    The `value` of this widget accepts a byte string.  The byte string is the
    raw video data that you want the browser to display.  You can explicitly
    define the format of the byte string using the `format` trait (which
    defaults to "mp4").

    If you pass `"url"` to the `"format"` trait, `value` will be interpreted
    as a URL as bytes encoded in UTF-8.
    """
    _view_name = Unicode('VideoView').tag(sync=True)
    _model_name = Unicode('VideoModel').tag(sync=True)

    # Define the custom state properties to sync with the front-end
    format = Unicode('mp4', help="The format of the video.").tag(sync=True)
    width = CUnicode(help="Width of the video in pixels.").tag(sync=True)
    height = CUnicode(help="Height of the video in pixels.").tag(sync=True)
    autoplay = Bool(True, help="When true, the video starts when it's displayed").tag(sync=True)
    loop = Bool(True, help="When true, the video will start from the beginning after finishing").tag(sync=True)
    controls = Bool(True, help="Specifies that video controls should be displayed (such as a play/pause button etc)").tag(sync=True)

    @classmethod
    def from_file(cls, filename, **kwargs):
        return cls._from_file('video', filename, **kwargs)

    def __repr__(self):
        return self._get_repr(Video)


@register
class Audio(_Media):
    """Displays a audio as a widget.

    The `value` of this widget accepts a byte string.  The byte string is the
    raw audio data that you want the browser to display.  You can explicitly
    define the format of the byte string using the `format` trait (which
    defaults to "mp3").

    If you pass `"url"` to the `"format"` trait, `value` will be interpreted
    as a URL as bytes encoded in UTF-8.
    """
    _view_name = Unicode('AudioView').tag(sync=True)
    _model_name = Unicode('AudioModel').tag(sync=True)

    # Define the custom state properties to sync with the front-end
    format = Unicode('mp3', help="The format of the audio.").tag(sync=True)
    autoplay = Bool(True, help="When true, the audio starts when it's displayed").tag(sync=True)
    loop = Bool(True, help="When true, the audio will start from the beginning after finishing").tag(sync=True)
    controls = Bool(True, help="Specifies that audio controls should be displayed (such as a play/pause button etc)").tag(sync=True)

    @classmethod
    def from_file(cls, filename, **kwargs):
        return cls._from_file('audio', filename, **kwargs)

    def __repr__(self):
        return self._get_repr(Audio)
