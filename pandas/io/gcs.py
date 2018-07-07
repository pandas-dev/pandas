""" GCS support for remote file interactivity """
try:
    import gcsfs
except ImportError:
    raise ImportError("The gcsfs library is required to handle GCS files")


def get_filepath_or_buffer(filepath_or_buffer, encoding=None,
                           compression=None, mode=None):

    if mode is None:
        mode = 'rb'

    fs = gcsfs.GCSFileSystem()
    filepath_or_buffer = fs.open(filepath_or_buffer, mode)
    return filepath_or_buffer, None, compression, True
