def GCSMap(root, gcs=None, check=False, create=False):
    """For backward compatibility"""
    import gcsfs

    gcs = gcs or gcsfs.GCSFileSystem.current()
    return gcs.get_mapper(root, check=check, create=create)
