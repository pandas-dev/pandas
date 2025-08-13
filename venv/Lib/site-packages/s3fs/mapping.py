from .core import S3FileSystem


def S3Map(root, s3, check=False, create=False):
    """Mirror previous class, not implemented in fsspec"""
    s3 = s3 or S3FileSystem.current()
    return s3.get_mapper(root, check=check, create=create)
