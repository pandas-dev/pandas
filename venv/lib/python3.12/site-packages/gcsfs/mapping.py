from .core import GCSFileSystem


def GCSMap(root, gcs=None, check=False, create=False):
    """For backward compatibility"""
    gcs = gcs or GCSFileSystem.current()
    return gcs.get_mapper(root, check=check, create=create)
